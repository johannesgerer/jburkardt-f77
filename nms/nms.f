      SUBROUTINE BAKSLD(NR,N,A,X,B)

c*********************************************************************72
C
C PURPOSE
C -------
C SOLVE  AX=B  WHERE A IS UPPER TRIANGULAR MATRIX.
C NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND
C THAT THIS ROUTINE TAKES ITS TRANSPOSE IMPLICITLY.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED)
C X(N)        <--  SOLUTION VECTOR
C B(N)         --> RIGHT-HAND SIDE VECTOR
C
C NOTE
C ----
C IF B IS NO LONGER REQUIRED BY CALLING ROUTINE,
C THEN VECTORS B AND X MAY SHARE THE SAME STORAGE.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(NR,1),X(N),B(N)
C
C SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE)
C
      I=N
      X(I)=B(I)/A(I,I)
      IF(N.EQ.1) RETURN
   30 IP1=I
      I=I-1
      SUM=0.D0
      DO 40 J=IP1,N
        SUM=SUM+A(J,I)*X(J)
   40 CONTINUE
      X(I)=(B(I)-SUM)/A(I,I)
      IF(I.GT.1) GO TO 30
      RETURN
      END
      SUBROUTINE CHLDCD(NR,N,A,DIAGMX,TOL,ADDMAX)

c*********************************************************************72
C
C PURPOSE
C -------
C FIND THE PERTURBED L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION
C OF A+D, WHERE D IS A NON-NEGATIVE DIAGONAL MATRIX ADDED TO A IF
C NECESSARY TO ALLOW THE CHOLESKY DECOMPOSITION TO CONTINUE.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--> ON ENTRY: MATRIX FOR WHICH TO FIND PERTURBED
C                       CHOLESKY DECOMPOSITION
C                  ON EXIT:  CONTAINS L OF LL+ DECOMPOSITION
C                  IN LOWER TRIANGULAR PART AND DIAGONAL OF "A"
C DIAGMX       --> MAXIMUM DIAGONAL ELEMENT OF "A"
C TOL          --> TOLERANCE
C ADDMAX      <--  MAXIMUM AMOUNT IMPLICITLY ADDED TO DIAGONAL OF "A"
C                  IN FORMING THE CHOLESKY DECOMPOSITION OF A+D
C INTERNAL VARIABLES
C ------------------
C AMINL    SMALLEST ELEMENT ALLOWED ON DIAGONAL OF L
C AMNLSQ   =AMINL**2
C OFFMAX   MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN OF A
C
C
C DESCRIPTION
C -----------
C THE NORMAL CHOLESKY DECOMPOSITION IS PERFORMED.  HOWEVER, IF AT ANY
C POINT THE ALGORITHM WOULD ATTEMPT TO SET L(I,I)=SQRT(TEMP)
C WITH TEMP < TOL*DIAGMX, THEN L(I,I) IS SET TO SQRT(TOL*DIAGMX)
C INSTEAD.  THIS IS EQUIVALENT TO ADDING TOL*DIAGMX-TEMP TO A(I,I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(NR,1)
C
      ADDMAX=0.D0
      AMINL=SQRT(DIAGMX*TOL)
      AMNLSQ=AMINL*AMINL
C
C FORM COLUMN J OF L
C
      DO 100 J=1,N
C FIND DIAGONAL ELEMENTS OF L
        SUM=0.D0
        IF(J.EQ.1) GO TO 20
        JM1=J-1
        DO 10 K=1,JM1
          SUM=SUM + A(J,K)*A(J,K)
   10   CONTINUE
   20   TEMP=A(J,J)-SUM
        IF(TEMP.LT.AMNLSQ) GO TO 30
C       IF(TEMP.GE.AMINL**2)
C       THEN
          A(J,J)=SQRT(TEMP)
          GO TO 40
C       ELSE
C
C FIND MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN
   30     OFFMAX=0.D0
          IF(J.EQ.N) GO TO 37
          JP1=J+1
          DO 35 I=JP1,N
            IF(ABS(A(I,J)).GT.OFFMAX) OFFMAX=ABS(A(I,J))
   35     CONTINUE
   37     IF(OFFMAX.LE.AMNLSQ) OFFMAX=AMNLSQ
C
C ADD TO DIAGONAL ELEMENT  TO ALLOW CHOLESKY DECOMPOSITION TO CONTINUE
          A(J,J)=SQRT(OFFMAX)
          ADDMAX=MAX(ADDMAX,OFFMAX-TEMP)
C       ENDIF
C
C FIND I,J ELEMENT OF LOWER TRIANGULAR MATRIX
   40   IF(J.EQ.N) GO TO 100
        JP1=J+1
        DO 70 I=JP1,N
          SUM=0.0D0
          IF(J.EQ.1) GO TO 60
          JM1=J-1
          DO 50 K=1,JM1
            SUM=SUM+A(I,K)*A(J,K)
   50     CONTINUE
   60     A(I,J)=(A(I,J)-SUM)/A(J,J)
   70   CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CHLHSD(NR,N,A,EPSM,SX,UDIAG)

c*********************************************************************72
c
C PURPOSE
C -------
C FIND THE L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION OF THE PERTURBED
C MODEL HESSIAN MATRIX A+MU*I(WHERE MU\0 AND I IS THE IDENTITY MATRIX)
C WHICH IS SAFELY POSITIVE DEFINITE.  IF A IS SAFELY POSITIVE DEFINITE
C UPON ENTRY, THEN MU=0.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--> ON ENTRY; "A" IS MODEL HESSIAN (ONLY LOWER
C                  TRIANGULAR PART AND DIAGONAL STORED)
C                  ON EXIT:  A CONTAINS L OF LL+ DECOMPOSITION OF
C                  PERTURBED MODEL HESSIAN IN LOWER TRIANGULAR
C                  PART AND DIAGONAL AND CONTAINS HESSIAN IN UPPER
C                  TRIANGULAR PART AND UDIAG
C EPSM         --> MACHINE EPSILON
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C UDIAG(N)    <--  ON EXIT: CONTAINS DIAGONAL OF HESSIAN
C
C INTERNAL VARIABLES
C ------------------
C TOL              TOLERANCE
C DIAGMN           MINIMUM ELEMENT ON DIAGONAL OF A
C DIAGMX           MAXIMUM ELEMENT ON DIAGONAL OF A
C OFFMAX           MAXIMUM OFF-DIAGONAL ELEMENT OF A
C OFFROW           SUM OF OFF-DIAGONAL ELEMENTS IN A ROW OF A
C EVMIN            MINIMUM EIGENVALUE OF A
C EVMAX            MAXIMUM EIGENVALUE OF A
C
C DESCRIPTION
C -----------
C 1. IF "A" HAS ANY NEGATIVE DIAGONAL ELEMENTS, THEN CHOOSE MU>0
C SUCH THAT THE DIAGONAL OF A:=A+MU*I IS ALL POSITIVE
C WITH THE RATIO OF ITS SMALLEST TO LARGEST ELEMENT ON THE
C ORDER OF SQRT(EPSM).
C
C 2. "A" UNDERGOES A PERTURBED CHOLESKY DECOMPOSITION WHICH
C RESULTS IN AN LL+ DECOMPOSITION OF A+D, WHERE D IS A
C NON-NEGATIVE DIAGONAL MATRIX WHICH IS IMPLICITLY ADDED TO
C "A" DURING THE DECOMPOSITION IF "A" IS NOT POSITIVE DEFINITE.
C "A" IS RETAINED AND NOT CHANGED DURING THIS PROCESS BY
C COPYING L INTO THE UPPER TRIANGULAR PART OF "A" AND THE
C DIAGONAL INTO UDIAG.  THEN THE CHOLESKY DECOMPOSITION ROUTINE
C IS CALLED.  ON RETURN, ADDMAX CONTAINS MAXIMUM ELEMENT OF D.
C
C 3. IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2
C AND RETURN IS MADE TO CALLING PROGRAM.  OTHERWISE,
C THE MINIMUM NUMBER SDD WHICH MUST BE ADDED TO THE
C DIAGONAL OF A TO MAKE IT SAFELY STRICTLY DIAGONALLY DOMINANT
C IS CALCULATED.  SINCE A+ADDMAX*I AND A+SDD*I ARE SAFELY
C POSITIVE DEFINITE, CHOOSE MU=MIN(ADDMAX,SDD) AND DECOMPOSE
C A+MU*I TO OBTAIN L.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION A(NR,1),SX(N),UDIAG(N)
C
C SCALE HESSIAN
C PRE- AND POST- MULTIPLY "A" BY INV(SX)
C
      DO 20 J=1,N
        DO 10 I=J,N
          A(I,J)=A(I,J)/(SX(I)*SX(J))
   10   CONTINUE
   20 CONTINUE
C
C STEP1
C -----
C NOTE:  IF A DIFFERENT TOLERANCE IS DESIRED THROUGHOUT THIS
C ALGORITHM, CHANGE TOLERANCE HERE:
      TOL=SQRT(EPSM)
C
      DIAGMX=A(1,1)
      DIAGMN=A(1,1)
      IF(N.EQ.1) GO TO 35
      DO 30 I=2,N
        IF(A(I,I).LT.DIAGMN) DIAGMN=A(I,I)
        IF(A(I,I).GT.DIAGMX) DIAGMX=A(I,I)
   30 CONTINUE
   35 POSMAX=MAX(DIAGMX,0.D0)
C
C DIAGMN .LE. 0
C
      IF(DIAGMN.GT.POSMAX*TOL) GO TO 100
C     IF(DIAGMN.LE.POSMAX*TOL)
C     THEN
        AMU=TOL*(POSMAX-DIAGMN)-DIAGMN
        IF(AMU.NE.0.D0) GO TO 60
C       IF(AMU.EQ.0.)
C       THEN
C
C FIND LARGEST OFF-DIAGONAL ELEMENT OF A
          OFFMAX=0.D0
          IF(N.EQ.1) GO TO 50
          DO 45 I=2,N
            IM1=I-1
            DO 40 J=1,IM1
              IF(ABS(A(I,J)).GT.OFFMAX) OFFMAX=ABS(A(I,J))
   40       CONTINUE
   45     CONTINUE
   50     AMU=OFFMAX
          IF(AMU.NE.0.D0) GO TO 55
C         IF(AMU.EQ.0.)
C         THEN
            AMU=1.0D0
            GO TO 60
C         ELSE
   55       AMU=AMU*(1.0D0+TOL)
C         ENDIF
C       ENDIF
C
C A=A + MU*I
C
   60   DO 65 I=1,N
          A(I,I)=A(I,I)+AMU
   65   CONTINUE
        DIAGMX=DIAGMX+AMU
C     ENDIF
C
C STEP2
C -----
C COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART
C AND DIAGONAL OF "A" TO UDIAG
C
  100 CONTINUE
      DO 110 J=1,N
        UDIAG(J)=A(J,J)
        IF(J.EQ.N) GO TO 110
        JP1=J+1
        DO 105 I=JP1,N
          A(J,I)=A(I,J)
  105   CONTINUE
  110 CONTINUE
C
      CALL CHLDCD(NR,N,A,DIAGMX,TOL,ADDMAX)
C
C
C STEP3
C -----
C IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2,
C THE LL+ DECOMPOSITION HAS BEEN DONE, AND WE RETURN.
C OTHERWISE, ADDMAX>0.  PERTURB "A" SO THAT IT IS SAFELY
C DIAGONALLY DOMINANT AND FIND LL+ DECOMPOSITION
C
      IF(ADDMAX.LE.0.D0) GO TO 170
C     IF(ADDMAX.GT.0.)
C     THEN
C
C RESTORE ORIGINAL "A" (LOWER TRIANGULAR PART AND DIAGONAL)
C
        DO 120 J=1,N
          A(J,J)=UDIAG(J)
          IF(J.EQ.N) GO TO 120
          JP1=J+1
          DO 115 I=JP1,N
            A(I,J)=A(J,I)
  115     CONTINUE
  120   CONTINUE
C
C FIND SDD SUCH THAT A+SDD*I IS SAFELY POSITIVE DEFINITE
C NOTE:  EVMIN<0 SINCE A IS NOT POSITIVE DEFINITE;
C
        EVMIN=0.D0
        EVMAX=A(1,1)
        DO 150 I=1,N
          OFFROW=0.D0
          IF(I.EQ.1) GO TO 135
          IM1=I-1
          DO 130 J=1,IM1
            OFFROW=OFFROW+ABS(A(I,J))
  130     CONTINUE
  135     IF(I.EQ.N) GO TO 145
          IP1=I+1
          DO 140 J=IP1,N
            OFFROW=OFFROW+ABS(A(J,I))
  140     CONTINUE
  145     EVMIN=MIN(EVMIN,A(I,I)-OFFROW)
          EVMAX=MAX(EVMAX,A(I,I)+OFFROW)
  150   CONTINUE
        SDD=TOL*(EVMAX-EVMIN)-EVMIN
C
C PERTURB "A" AND DECOMPOSE AGAIN
C
        AMU=MIN(SDD,ADDMAX)
        DO 160 I=1,N
          A(I,I)=A(I,I)+AMU
          UDIAG(I)=A(I,I)
  160   CONTINUE
C
C "A" NOW GUARANTEED SAFELY POSITIVE DEFINITE
C
        CALL CHLDCD(NR,N,A,0.0D0,TOL,ADDMAX)
C     ENDIF
C
C UNSCALE HESSIAN AND CHOLESKY DECOMPOSITION MATRIX
C
  170 DO 190 J=1,N
        DO 175 I=J,N
          A(I,J)=SX(I)*A(I,J)
  175   CONTINUE
        IF(J.EQ.1) GO TO 185
        JM1=J-1
        DO 180 I=1,JM1
          A(I,J)=SX(I)*SX(J)*A(I,J)
  180   CONTINUE
  185   UDIAG(J)=UDIAG(J)*SX(J)*SX(J)
  190 CONTINUE
      RETURN
      END
      SUBROUTINE D1FCND(N,X,G)

c*********************************************************************72
c
C PURPOSE
C -------
C DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC
C WHEN SPECIFIC ANALYTIC GRADIENT FUNCTION NOT SUPPLIED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION X(N),G(N)
      G(N)=G(N)
      X(N)=X(N)
      STOP
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)

c*********************************************************************72
c
C***BEGIN PROLOGUE  D1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  831014   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns double precision machine dependent constants
C***DESCRIPTION
C     From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C
C     D1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subprogram with one (input) argument, and can be called
C     as follows, for example
C
C          D = D1MACH(I)
C
C     where I=1,...,5.  The (output) value of D above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  Double-precision machine constants
C  D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C  D1MACH( 3) = B**(-T), the smallest relative spacing.
C  D1MACH( 4) = B**(1-T), the largest relative spacing.
C  D1MACH( 5) = LOG10(B)
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
C
C      DATA SMALL(1) / O"00604000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C
C      DATA LARGE(1) / O"37767777777777777777" /
C      DATA LARGE(2) / O"37167777777777777777" /
C
C      DATA RIGHT(1) / O"15604000000000000000" /
C      DATA RIGHT(2) / O"15000000000000000000" /
C
C      DATA DIVER(1) / O"15614000000000000000" /
C      DATA DIVER(2) / O"15010000000000000000" /
C
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
C
C     DATA SMALL(1) / X'9000400000000000' /
C     DATA SMALL(2) / X'8FD1000000000000' /
C
C     DATA LARGE(1) / X'6FFF7FFFFFFFFFFF' /
C     DATA LARGE(2) / X'6FD07FFFFFFFFFFF' /
C
C     DATA RIGHT(1) / X'FF74400000000000' /
C     DATA RIGHT(2) / X'FF45000000000000' /
C
C     DATA DIVER(1) / X'FF75400000000000' /
C     DATA DIVER(2) / X'FF46000000000000' /
C
C     DATA LOG10(1) / X'FFD04D104D427DE7' /
C     DATA LOG10(2) / X'FFA17DE623E2566A' /
C
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTATNS FOR THE IBM PC FAMILY (D. KAHANER NBS)
C
      DATA DMACH/2.23D-308,1.79D+308,1.11D-16,2.22D-16,
     *  0.301029995663981195D0/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1),LOG10(2) / "177464202324, "476747767461 /
C
C
C     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
C       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
C       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
C
C      DATA SMALL(1),SMALL(2) / X'00100000', X'00000000' /
C      DATA LARGE(1),LARGE(2) / X'7FEFFFFF', X'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / X'3CA00000', X'00000000' /
C      DATA DIVER(1),DIVER(2) / X'3CB00000', X'00000000' /
C      DATA LOG10(1),LOG10(2) / X'3FD34413', X'509F79FF' /
C
C
C     MACHINE CONSTANTS FOR VAX 11/780
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C   MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING)
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C    *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C    ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS***
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR( 'D1MACH -- I OUT OF BOUNDS',25,1,2)

      D1MACH = DMACH(I)

      RETURN
      END
      SUBROUTINE D1MPYQ(M,N,A,LDA,V,W)

c*********************************************************************72
c
C***BEGIN PROLOGUE  D1MPYQ
C***REFER TO  DNSQ,DNSQE
C
C     SUBROUTINE D1MPYQ
C
C     Given an M by N matrix A, this subroutine computes A*Q where
C     Q is the product of 2*(N - 1) transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
C     eliminate elements in the I-th and N-th planes, respectively.
C     Q itself is not given, rather the information to recover the
C     GV, GW rotations is supplied.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE D1MPYQ(M,N,A,LDA,V,W)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A.
C
C       N IS a positive integer input variable set to the number
C         of columns of A.
C
C       A is an M by N array. On input A must contain the matrix
C         to be postmultiplied by the orthogonal matrix Q
C         described above. On output A*Q has replaced A.
C
C       LDA is a positive integer input variable not less than M
C         which specifies the leading dimension of the array A.
C
C       V is an input array of length N. V(I) must contain the
C         information necessary to recover the Givens rotation GV(I)
C         described above.
C
C       W is an input array of length N. W(I) must contain the
C         information necessary to recover the Givens rotation GW(I)
C         described above.
C
C     Subroutines called
C
C       FORTRAN-supplied ... DABS,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  D1MPYQ
      DOUBLE PRECISION DABS, DSQRT
      INTEGER I, J, LDA, M, N, NM1, NMJ
      DOUBLE PRECISION A(LDA,N), COS, ONE, SIN, TEMP, V(N), W(N)
      SAVE ONE
      DATA ONE /1.0D0/
C
C     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
C
C***FIRST EXECUTABLE STATEMENT  D1MPYQ
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (DABS(V(J)) .GT. ONE) COS = ONE/V(J)
         IF (DABS(V(J)) .GT. ONE) SIN = DSQRT(ONE-COS**2)
         IF (DABS(V(J)) .LE. ONE) SIN = V(J)
         IF (DABS(V(J)) .LE. ONE) COS = DSQRT(ONE-SIN**2)
         DO 10 I = 1, M
            TEMP = COS*A(I,J) - SIN*A(I,N)
            A(I,N) = SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
C
      DO 40 J = 1, NM1
         IF (DABS(W(J)) .GT. ONE) COS = ONE/W(J)
         IF (DABS(W(J)) .GT. ONE) SIN = DSQRT(ONE-COS**2)
         IF (DABS(W(J)) .LE. ONE) SIN = W(J)
         IF (DABS(W(J)) .LE. ONE) COS = DSQRT(ONE-SIN**2)
         DO 30 I = 1, M
            TEMP = COS*A(I,J) + SIN*A(I,N)
            A(I,N) = -SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)

c*********************************************************************72
c
C***BEGIN PROLOGUE  D1UPDT
C***REFER TO  DNSQ,DNSQE
C
C     SUBROUTINE D1UPDT
C
C     Given an M by N lower trapezoidal matrix S, an M-vector U,
C     and an N-vector V, the problem is to determine an
C     orthogonal matrix Q such that
C
C                   t
C           (S + U*V )*Q
C
C     is again lower trapezoidal.
C
C     This subroutine determines Q as the product of 2*(N - 1)
C     transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     where GV(I), GW(I) are Givens rotations in the (I,N) plane
C     which eliminate elements in the I-th and N-th planes,
C     respectively. Q itself is not accumulated, rather the
C     information to recover the GV, GW rotations is returned.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of S.
C
C       N is a positive integer input variable set to the number
C         of columns of S. N must not exceed M.
C
C       S is an array of length LS. On input S must contain the lower
C         trapezoidal matrix S stored by columns. On output S contains
C         the lower trapezoidal matrix produced as described above.
C
C       LS is a positive integer input variable not less than
C         (N*(2*M-N+1))/2.
C
C       U is an input array of length M which must contain the
C         vector U.
C
C       V is an array of length N. On input V must contain the vector
C         V. On output V(I) contains the information necessary to
C         recover the Givens rotation GV(I) described above.
C
C       W is an output array of length M. W(I) contains information
C         necessary to recover the Givens rotation GW(I) described
C         above.
C
C       SING is a LOGICAL output variable. SING is set TRUE if any
C         of the diagonal elements of the output S are zero. Otherwise
C         SING is set FALSE.
C
C     Subprograms called
C
C       SLATEC-supplied ... D1MACH
C
C       FORTRAN-supplied ... DABS,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More,
C     John L. Nazareth
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  D1UPDT
      DOUBLE PRECISION D1MACH
      DOUBLE PRECISION DABS, DSQRT
      INTEGER I, J, JJ, L, LS, M, N, NM1, NMJ
      DOUBLE PRECISION COS, COTAN, GIANT, ONE, P25, P5, S(LS),
     1     SIN, TAN, TAU, TEMP, U(M), V(N), W(M), ZERO
      LOGICAL SING
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
C
C     GIANT IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  D1UPDT
      GIANT = D1MACH(2)
C
C     INITIALIZE THE DIAGONAL ELEMENT POINTER.
C
      JJ = (N*(2*M - N + 1))/2 - (M - N)
C
C     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
C
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
C
C     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
C     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (V(J) .EQ. ZERO) GO TO 50
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF V.
C
         IF (DABS(V(N)) .GE. DABS(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            SIN = P5/DSQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (DABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 30
   20    CONTINUE
            TAN = V(J)/V(N)
            COS = P5/DSQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
   30    CONTINUE
C
C        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
C        NECESSARY TO RECOVER THE GIVENS ROTATION.
C
         V(N) = SIN*V(J) + COS*V(N)
         V(J) = TAU
C
C        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
C
         L = JJ
         DO 40 I = J, M
            TEMP = COS*S(L) - SIN*W(I)
            W(I) = SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
C
C     ELIMINATE THE SPIKE.
C
      SING = .FALSE.
      IF (NM1 .LT. 1) GO TO 140
      DO 130 J = 1, NM1
         IF (W(J) .EQ. ZERO) GO TO 120
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF THE SPIKE.
C
         IF (DABS(S(JJ)) .GE. DABS(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            SIN = P5/DSQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (DABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 100
   90    CONTINUE
            TAN = W(J)/S(JJ)
            COS = P5/DSQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
  100    CONTINUE
C
C        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
C
         L = JJ
         DO 110 I = J, M
            TEMP = COS*S(L) + SIN*W(I)
            W(I) = -SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
C
C        STORE THE INFORMATION NECESSARY TO RECOVER THE
C        GIVENS ROTATION.
C
         W(J) = TAU
  120    CONTINUE
C
C        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
C
         IF (S(JJ) .EQ. ZERO) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
C
C     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
C
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (S(JJ) .EQ. ZERO) SING = .TRUE.
      RETURN
      END
      SUBROUTINE D2FCND(NR,N,X,H)

c*********************************************************************72
c
C PURPOSE
C -------
C DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC
C WHEN SPECIFIC ANALYTIC HESSIAN FUNCTION NOT SUPPLIED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION X(N),H(NR,1)
      H(NR,1)=H(NR,1)
      X(N)=X(N)
      STOP
      END
      DOUBLE PRECISION FUNCTION D9LGMC(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  D9LGMC
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C7E
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(R9LGMC-S D9LGMC-D C9LGMC-C),
C             COMPLETE GAMMA FUNCTION,CORRECTION TERM,GAMMA FUNCTION,
C             LOG GAMMA,LOGARITHM,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the  d.p. log Gamma correction factor for
C            X .GE. 10. so that DLOG(DGAMMA(X)) = DLOG(DSQRT(2*PI)) +
C            (X-5.)*DLOG(X) - X + D9LGMC(X)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C DLOG (DGAMMA(X)) = DLOG(DSQRT(2*PI)) + (X-.5)*DLOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS,XERROR
C***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      SAVE ALGMCS, NALGM, XBIG, XMAX
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
      DATA NALGM, XBIG, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (NALGM.NE.0) GO TO 10
      NALGM = INITDS (ALGMCS, 15, SNGL(D1MACH(3)) )
      XBIG = 1.0D0/DSQRT(D1MACH(3))
      XMAX = DEXP (DMIN1(DLOG(D1MACH(2)/12.D0), -DLOG(12.D0*D1MACH(1))))
C
 10   IF (X.LT.10.D0) CALL XERROR ( 'D9LGMC  X MUST BE GE 10', 23, 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
C
 20   D9LGMC = 0.D0
      CALL XERROR ( 'D9LGMC  X SO BIG D9LGMC UNDERFLOWS', 34, 2, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DASUM
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A3A
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SASUM-S DASUM-D SCASUM-C),ADD,
C             LINEAR ALGEBRA,MAGNITUDE,SUM,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Sum of magnitudes of d.p. vector components
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DASUM  double precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of double precision DX.
C     DASUM = sum from 0 to N-1 of DABS(DX(1+I*INCX))
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DASUM
C
      DOUBLE PRECISION DX(1)
C***FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          DASUM = DASUM + DABS(DX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DASUM = DASUM + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
         DASUM = DASUM + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2))
     1   + DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DASYJY(FUNJY,X,FNU,FLGJY,IN,Y,WK,IFLW)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DASYJY
C***REFER TO  DBESJ,DBESY
C***ROUTINES CALLED  D1MACH,I1MACH
C***DESCRIPTION
C
C                 DASYJY computes Bessel functions J and Y
C               for arguments X.GT.0.0 and orders FNU .GE. 35.0
C               on FLGJY = 1 and FLGJY = -1 respectively
C
C                                  INPUT
C
C      FUNJY - External function JAIRY or YAIRY
C          X - Argument, X.GT.0.0D0
C        FNU - Order of the first Bessel function
C      FLGJY - Selection flag
C              FLGJY =  1.0D0 gives the J function
C              FLGJY = -1.0D0 gives the Y function
C         IN - Number of functions desired, IN = 1 or 2
C
C                                  OUTPUT
C
C         Y  - A vector whose first IN components contain the sequence
C       IFLW - A flag indicating underflow or overflow
C                    return variables for BESJ only
C      WK(1) = 1 - (X/FNU)**2 = W**2
C      WK(2) = DSQRT(DABS(WK(1)))
C      WK(3) = DABS(WK(2) - DATAN(WK(2)))  or
C              DABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
C            = DABS((2/3)*ZETA**(3/2))
C      WK(4) = FNU*WK(3)
C      WK(5) = (1.5*WK(3)*FNU)**(1/3) = DSQRT(ZETA)*FNU**(1/3)
C      WK(6) = DSIGN(1.,W**2)*WK(5)**2 = DSIGN(1.,W**2)*ZETA*FNU**(2/3)
C      WK(7) = FNU**(1/3)
C
C                                  Written by
C                                  D. E. Amos
C
C     Abstract   **** A Double Precision Routine ****
C         DASYJY implements the uniform asymptotic expansion of
C         the J and Y Bessel functions for FNU.GE.35 and real
C         X.GT.0.0D0. The forms are identical except for a change
C         in sign of some of the terms. This change in sign is
C         accomplished by means of the flag FLGJY = 1 or -1. On
C         FLGJY = 1 the Airy functions AI(X) and DAI(X) are
C         supplied by the external function JAIRY, and on
C         FLGJY = -1 the Airy functions BI(X) and DBI(X) are
C         supplied by the external funtion YAIRY.
C***END PROLOGUE  DASYJY
      INTEGER I, IFLW, IN, J, JN,JR,JU,K, KB,KLAST,KMAX,KP1, KS, KSP1,
     * KSTEMP, L, LR, LRP1, ISETA, ISETB
      INTEGER I1MACH
      DOUBLE PRECISION ABW2, AKM, ALFA, ALFA1, ALFA2, AP, AR, ASUM, AZ,
     * BETA, BETA1, BETA2, BETA3, BR, BSUM, C, CON1, CON2,
     * CON548,CR,CRZ32, DFI,ELIM, DR,FI, FLGJY, FN, FNU,
     * FN2, GAMA, PHI,  RCZ, RDEN, RELB, RFN2,  RTZ, RZDEN,
     * SA, SB, SUMA, SUMB, S1, TA, TAU, TB, TFN, TOL, TOLS, T2, UPOL,
     *  WK, X, XX, Y, Z, Z32
      DOUBLE PRECISION D1MACH
      DIMENSION Y(*), WK(*), C(65)
      DIMENSION ALFA(26,4), BETA(26,5)
      DIMENSION ALFA1(26,2), ALFA2(26,2)
      DIMENSION BETA1(26,2), BETA2(26,2), BETA3(26,1)
      DIMENSION GAMA(26), KMAX(5), AR(8), BR(10), UPOL(10)
      DIMENSION CR(10), DR(10)
      EQUIVALENCE (ALFA(1,1),ALFA1(1,1))
      EQUIVALENCE (ALFA(1,3),ALFA2(1,1))
      EQUIVALENCE (BETA(1,1),BETA1(1,1))
      EQUIVALENCE (BETA(1,3),BETA2(1,1))
      EQUIVALENCE (BETA(1,5),BETA3(1,1))
      SAVE TOLS,CON1, CON2, CON548, AR, BR, C,
     1 ALFA1, ALFA2, BETA1, BETA2, BETA3, GAMA
      DATA TOLS            /-6.90775527898214D+00/
      DATA CON1,CON2,CON548/
     1 6.66666666666667D-01, 3.33333333333333D-01, 1.04166666666667D-01/
      DATA  AR(1),  AR(2),  AR(3),  AR(4),  AR(5),  AR(6),  AR(7),
     A      AR(8)          / 8.35503472222222D-02, 1.28226574556327D-01,
     1 2.91849026464140D-01, 8.81627267443758D-01, 3.32140828186277D+00,
     2 1.49957629868626D+01, 7.89230130115865D+01, 4.74451538868264D+02/
      DATA  BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),
     A      BR(9), BR(10)  /-1.45833333333333D-01,-9.87413194444444D-02,
     1-1.43312053915895D-01,-3.17227202678414D-01,-9.42429147957120D-01,
     2-3.51120304082635D+00,-1.57272636203680D+01,-8.22814390971859D+01,
     3-4.92355370523671D+02,-3.31621856854797D+03/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),
     1     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),
     2     C(19), C(20), C(21), C(22), C(23), C(24)/
     3       -2.08333333333333D-01,        1.25000000000000D-01,
     4        3.34201388888889D-01,       -4.01041666666667D-01,
     5        7.03125000000000D-02,       -1.02581259645062D+00,
     6        1.84646267361111D+00,       -8.91210937500000D-01,
     7        7.32421875000000D-02,        4.66958442342625D+00,
     8       -1.12070026162230D+01,        8.78912353515625D+00,
     9       -2.36408691406250D+00,        1.12152099609375D-01,
     A       -2.82120725582002D+01,        8.46362176746007D+01,
     B       -9.18182415432400D+01,        4.25349987453885D+01,
     C       -7.36879435947963D+00,        2.27108001708984D-01,
     D        2.12570130039217D+02,       -7.65252468141182D+02,
     E        1.05999045252800D+03,       -6.99579627376133D+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),
     1     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),
     2     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/
     3        2.18190511744212D+02,       -2.64914304869516D+01,
     4        5.72501420974731D-01,       -1.91945766231841D+03,
     5        8.06172218173731D+03,       -1.35865500064341D+04,
     6        1.16553933368645D+04,       -5.30564697861340D+03,
     7        1.20090291321635D+03,       -1.08090919788395D+02,
     8        1.72772750258446D+00,        2.02042913309661D+04,
     9       -9.69805983886375D+04,        1.92547001232532D+05,
     A       -2.03400177280416D+05,        1.22200464983017D+05,
     B       -4.11926549688976D+04,        7.10951430248936D+03,
     C       -4.93915304773088D+02,        6.07404200127348D+00,
     D       -2.42919187900551D+05,        1.31176361466298D+06,
     E       -2.99801591853811D+06,        3.76327129765640D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),
     1     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),
     2     C(65)/
     3       -2.81356322658653D+06,        1.26836527332162D+06,
     4       -3.31645172484564D+05,        4.52187689813627D+04,
     5       -2.49983048181121D+03,        2.43805296995561D+01,
     6        3.28446985307204D+06,       -1.97068191184322D+07,
     7        5.09526024926646D+07,       -7.41051482115327D+07,
     8        6.63445122747290D+07,       -3.75671766607634D+07,
     9        1.32887671664218D+07,       -2.78561812808645D+06,
     A        3.08186404612662D+05,       -1.38860897537170D+04,
     B        1.10017140269247D+02/
      DATA ALFA1(1,1), ALFA1(2,1), ALFA1(3,1), ALFA1(4,1), ALFA1(5,1),
     1     ALFA1(6,1), ALFA1(7,1), ALFA1(8,1), ALFA1(9,1), ALFA1(10,1),
     2     ALFA1(11,1),ALFA1(12,1),ALFA1(13,1),ALFA1(14,1),ALFA1(15,1),
     3     ALFA1(16,1),ALFA1(17,1),ALFA1(18,1),ALFA1(19,1),ALFA1(20,1),
     4     ALFA1(21,1),ALFA1(22,1),ALFA1(23,1),ALFA1(24,1),ALFA1(25,1),
     5     ALFA1(26,1)     /-4.44444444444444D-03,-9.22077922077922D-04,
     6-8.84892884892885D-05, 1.65927687832450D-04, 2.46691372741793D-04,
     7 2.65995589346255D-04, 2.61824297061501D-04, 2.48730437344656D-04,
     8 2.32721040083232D-04, 2.16362485712365D-04, 2.00738858762752D-04,
     9 1.86267636637545D-04, 1.73060775917876D-04, 1.61091705929016D-04,
     1 1.50274774160908D-04, 1.40503497391270D-04, 1.31668816545923D-04,
     2 1.23667445598253D-04, 1.16405271474738D-04, 1.09798298372713D-04,
     3 1.03772410422993D-04, 9.82626078369363D-05, 9.32120517249503D-05,
     4 8.85710852478712D-05, 8.42963105715700D-05, 8.03497548407791D-05/
      DATA ALFA1(1,2), ALFA1(2,2), ALFA1(3,2), ALFA1(4,2), ALFA1(5,2),
     1     ALFA1(6,2), ALFA1(7,2), ALFA1(8,2), ALFA1(9,2), ALFA1(10,2),
     2     ALFA1(11,2),ALFA1(12,2),ALFA1(13,2),ALFA1(14,2),ALFA1(15,2),
     3     ALFA1(16,2),ALFA1(17,2),ALFA1(18,2),ALFA1(19,2),ALFA1(20,2),
     4     ALFA1(21,2),ALFA1(22,2),ALFA1(23,2),ALFA1(24,2),ALFA1(25,2),
     5     ALFA1(26,2)     / 6.93735541354589D-04, 2.32241745182922D-04,
     6-1.41986273556691D-05,-1.16444931672049D-04,-1.50803558053049D-04,
     7-1.55121924918096D-04,-1.46809756646466D-04,-1.33815503867491D-04,
     8-1.19744975684254D-04,-1.06184319207974D-04,-9.37699549891194D-05,
     9-8.26923045588193D-05,-7.29374348155221D-05,-6.44042357721016D-05,
     1-5.69611566009369D-05,-5.04731044303562D-05,-4.48134868008883D-05,
     2-3.98688727717599D-05,-3.55400532972042D-05,-3.17414256609022D-05,
     3-2.83996793904175D-05,-2.54522720634871D-05,-2.28459297164725D-05,
     4-2.05352753106481D-05,-1.84816217627666D-05,-1.66519330021394D-05/
      DATA ALFA2(1,1), ALFA2(2,1), ALFA2(3,1), ALFA2(4,1), ALFA2(5,1),
     1     ALFA2(6,1), ALFA2(7,1), ALFA2(8,1), ALFA2(9,1), ALFA2(10,1),
     2     ALFA2(11,1),ALFA2(12,1),ALFA2(13,1),ALFA2(14,1),ALFA2(15,1),
     3     ALFA2(16,1),ALFA2(17,1),ALFA2(18,1),ALFA2(19,1),ALFA2(20,1),
     4     ALFA2(21,1),ALFA2(22,1),ALFA2(23,1),ALFA2(24,1),ALFA2(25,1),
     5     ALFA2(26,1)     /-3.54211971457744D-04,-1.56161263945159D-04,
     6 3.04465503594936D-05, 1.30198655773243D-04, 1.67471106699712D-04,
     7 1.70222587683593D-04, 1.56501427608595D-04, 1.36339170977445D-04,
     8 1.14886692029825D-04, 9.45869093034688D-05, 7.64498419250898D-05,
     9 6.07570334965197D-05, 4.74394299290509D-05, 3.62757512005344D-05,
     1 2.69939714979225D-05, 1.93210938247939D-05, 1.30056674793963D-05,
     2 7.82620866744497D-06, 3.59257485819352D-06, 1.44040049814252D-07,
     3-2.65396769697939D-06,-4.91346867098486D-06,-6.72739296091248D-06,
     4-8.17269379678658D-06,-9.31304715093561D-06,-1.02011418798016D-05/
      DATA ALFA2(1,2), ALFA2(2,2), ALFA2(3,2), ALFA2(4,2), ALFA2(5,2),
     1     ALFA2(6,2), ALFA2(7,2), ALFA2(8,2), ALFA2(9,2), ALFA2(10,2),
     2     ALFA2(11,2),ALFA2(12,2),ALFA2(13,2),ALFA2(14,2),ALFA2(15,2),
     3     ALFA2(16,2),ALFA2(17,2),ALFA2(18,2),ALFA2(19,2),ALFA2(20,2),
     4     ALFA2(21,2),ALFA2(22,2),ALFA2(23,2),ALFA2(24,2),ALFA2(25,2),
     5     ALFA2(26,2)     / 3.78194199201773D-04, 2.02471952761816D-04,
     6-6.37938506318862D-05,-2.38598230603006D-04,-3.10916256027362D-04,
     7-3.13680115247576D-04,-2.78950273791323D-04,-2.28564082619141D-04,
     8-1.75245280340847D-04,-1.25544063060690D-04,-8.22982872820208D-05,
     9-4.62860730588116D-05,-1.72334302366962D-05, 5.60690482304602D-06,
     1 2.31395443148287D-05, 3.62642745856794D-05, 4.58006124490189D-05,
     2 5.24595294959114D-05, 5.68396208545815D-05, 5.94349820393104D-05,
     3 6.06478527578422D-05, 6.08023907788436D-05, 6.01577894539460D-05,
     4 5.89199657344698D-05, 5.72515823777593D-05, 5.52804375585853D-05/
      DATA BETA1(1,1), BETA1(2,1), BETA1(3,1), BETA1(4,1), BETA1(5,1),
     1     BETA1(6,1), BETA1(7,1), BETA1(8,1), BETA1(9,1), BETA1(10,1),
     2     BETA1(11,1),BETA1(12,1),BETA1(13,1),BETA1(14,1),BETA1(15,1),
     3     BETA1(16,1),BETA1(17,1),BETA1(18,1),BETA1(19,1),BETA1(20,1),
     4     BETA1(21,1),BETA1(22,1),BETA1(23,1),BETA1(24,1),BETA1(25,1),
     5     BETA1(26,1)     / 1.79988721413553D-02, 5.59964911064388D-03,
     6 2.88501402231133D-03, 1.80096606761054D-03, 1.24753110589199D-03,
     7 9.22878876572938D-04, 7.14430421727287D-04, 5.71787281789705D-04,
     8 4.69431007606482D-04, 3.93232835462917D-04, 3.34818889318298D-04,
     9 2.88952148495752D-04, 2.52211615549573D-04, 2.22280580798883D-04,
     1 1.97541838033063D-04, 1.76836855019718D-04, 1.59316899661821D-04,
     2 1.44347930197334D-04, 1.31448068119965D-04, 1.20245444949303D-04,
     3 1.10449144504599D-04, 1.01828770740567D-04, 9.41998224204238D-05,
     4 8.74130545753834D-05, 8.13466262162801D-05, 7.59002269646219D-05/
      DATA BETA1(1,2), BETA1(2,2), BETA1(3,2), BETA1(4,2), BETA1(5,2),
     1     BETA1(6,2), BETA1(7,2), BETA1(8,2), BETA1(9,2), BETA1(10,2),
     2     BETA1(11,2),BETA1(12,2),BETA1(13,2),BETA1(14,2),BETA1(15,2),
     3     BETA1(16,2),BETA1(17,2),BETA1(18,2),BETA1(19,2),BETA1(20,2),
     4     BETA1(21,2),BETA1(22,2),BETA1(23,2),BETA1(24,2),BETA1(25,2),
     5     BETA1(26,2)     /-1.49282953213429D-03,-8.78204709546389D-04,
     6-5.02916549572035D-04,-2.94822138512746D-04,-1.75463996970783D-04,
     7-1.04008550460816D-04,-5.96141953046458D-05,-3.12038929076098D-05,
     8-1.26089735980230D-05,-2.42892608575730D-07, 8.05996165414274D-06,
     9 1.36507009262147D-05, 1.73964125472926D-05, 1.98672978842134D-05,
     1 2.14463263790823D-05, 2.23954659232457D-05, 2.28967783814713D-05,
     2 2.30785389811178D-05, 2.30321976080909D-05, 2.28236073720349D-05,
     3 2.25005881105292D-05, 2.20981015361991D-05, 2.16418427448104D-05,
     4 2.11507649256221D-05, 2.06388749782171D-05, 2.01165241997082D-05/
      DATA BETA2(1,1), BETA2(2,1), BETA2(3,1), BETA2(4,1), BETA2(5,1),
     1     BETA2(6,1), BETA2(7,1), BETA2(8,1), BETA2(9,1), BETA2(10,1),
     2     BETA2(11,1),BETA2(12,1),BETA2(13,1),BETA2(14,1),BETA2(15,1),
     3     BETA2(16,1),BETA2(17,1),BETA2(18,1),BETA2(19,1),BETA2(20,1),
     4     BETA2(21,1),BETA2(22,1),BETA2(23,1),BETA2(24,1),BETA2(25,1),
     5     BETA2(26,1)     / 5.52213076721293D-04, 4.47932581552385D-04,
     6 2.79520653992021D-04, 1.52468156198447D-04, 6.93271105657044D-05,
     7 1.76258683069991D-05,-1.35744996343269D-05,-3.17972413350427D-05,
     8-4.18861861696693D-05,-4.69004889379141D-05,-4.87665447413787D-05,
     9-4.87010031186735D-05,-4.74755620890087D-05,-4.55813058138628D-05,
     1-4.33309644511266D-05,-4.09230193157750D-05,-3.84822638603221D-05,
     2-3.60857167535411D-05,-3.37793306123367D-05,-3.15888560772110D-05,
     3-2.95269561750807D-05,-2.75978914828336D-05,-2.58006174666884D-05,
     4-2.41308356761280D-05,-2.25823509518346D-05,-2.11479656768913D-05/
      DATA BETA2(1,2), BETA2(2,2), BETA2(3,2), BETA2(4,2), BETA2(5,2),
     1     BETA2(6,2), BETA2(7,2), BETA2(8,2), BETA2(9,2), BETA2(10,2),
     2     BETA2(11,2),BETA2(12,2),BETA2(13,2),BETA2(14,2),BETA2(15,2),
     3     BETA2(16,2),BETA2(17,2),BETA2(18,2),BETA2(19,2),BETA2(20,2),
     4     BETA2(21,2),BETA2(22,2),BETA2(23,2),BETA2(24,2),BETA2(25,2),
     5     BETA2(26,2)     /-4.74617796559960D-04,-4.77864567147321D-04,
     6-3.20390228067038D-04,-1.61105016119962D-04,-4.25778101285435D-05,
     7 3.44571294294968D-05, 7.97092684075675D-05, 1.03138236708272D-04,
     8 1.12466775262204D-04, 1.13103642108481D-04, 1.08651634848774D-04,
     9 1.01437951597662D-04, 9.29298396593364D-05, 8.40293133016090D-05,
     1 7.52727991349134D-05, 6.69632521975731D-05, 5.92564547323195D-05,
     2 5.22169308826976D-05, 4.58539485165361D-05, 4.01445513891487D-05,
     3 3.50481730031328D-05, 3.05157995034347D-05, 2.64956119950516D-05,
     4 2.29363633690998D-05, 1.97893056664022D-05, 1.70091984636413D-05/
      DATA BETA3(1,1), BETA3(2,1), BETA3(3,1), BETA3(4,1), BETA3(5,1),
     1     BETA3(6,1), BETA3(7,1), BETA3(8,1), BETA3(9,1), BETA3(10,1),
     2     BETA3(11,1),BETA3(12,1),BETA3(13,1),BETA3(14,1),BETA3(15,1),
     3     BETA3(16,1),BETA3(17,1),BETA3(18,1),BETA3(19,1),BETA3(20,1),
     4     BETA3(21,1),BETA3(22,1),BETA3(23,1),BETA3(24,1),BETA3(25,1),
     5     BETA3(26,1)     / 7.36465810572578D-04, 8.72790805146194D-04,
     6 6.22614862573135D-04, 2.85998154194304D-04, 3.84737672879366D-06,
     7-1.87906003636972D-04,-2.97603646594555D-04,-3.45998126832656D-04,
     8-3.53382470916038D-04,-3.35715635775049D-04,-3.04321124789040D-04,
     9-2.66722723047613D-04,-2.27654214122820D-04,-1.89922611854562D-04,
     1-1.55058918599094D-04,-1.23778240761874D-04,-9.62926147717644D-05,
     2-7.25178327714425D-05,-5.22070028895634D-05,-3.50347750511901D-05,
     3-2.06489761035552D-05,-8.70106096849767D-06, 1.13698686675100D-06,
     4 9.16426474122779D-06, 1.56477785428873D-05, 2.08223629482467D-05/
      DATA GAMA(1),   GAMA(2),   GAMA(3),   GAMA(4),   GAMA(5),
     1     GAMA(6),   GAMA(7),   GAMA(8),   GAMA(9),   GAMA(10),
     2     GAMA(11),  GAMA(12),  GAMA(13),  GAMA(14),  GAMA(15),
     3     GAMA(16),  GAMA(17),  GAMA(18),  GAMA(19),  GAMA(20),
     4     GAMA(21),  GAMA(22),  GAMA(23),  GAMA(24),  GAMA(25),
     5     GAMA(26)        / 6.29960524947437D-01, 2.51984209978975D-01,
     6 1.54790300415656D-01, 1.10713062416159D-01, 8.57309395527395D-02,
     7 6.97161316958684D-02, 5.86085671893714D-02, 5.04698873536311D-02,
     8 4.42600580689155D-02, 3.93720661543510D-02, 3.54283195924455D-02,
     9 3.21818857502098D-02, 2.94646240791158D-02, 2.71581677112934D-02,
     1 2.51768272973862D-02, 2.34570755306079D-02, 2.19508390134907D-02,
     2 2.06210828235646D-02, 1.94388240897881D-02, 1.83810633800683D-02,
     3 1.74293213231963D-02, 1.65685837786612D-02, 1.57865285987918D-02,
     4 1.50729501494096D-02, 1.44193250839955D-02, 1.38184805735342D-02/
C***FIRST EXECUTABLE STATEMENT  DASYJY
      TA = D1MACH(3)
      TOL = DMAX1(TA,1.0D-15)
      TB = D1MACH(5)
      JU = I1MACH(15)
      IF(FLGJY.EQ.1.0D0) GO TO 6
      JR = I1MACH(14)
      ELIM = 2.303D0*TB*(DBLE(FLOAT(-JU))-DBLE(FLOAT(JR)))
      GO TO 7
    6 CONTINUE
      ELIM = 2.303D0*(TB*DBLE(FLOAT(-JU))-3.0D0)
    7 CONTINUE
      FN = FNU
      IFLW = 0
      DO 170 JN=1,IN
        XX = X/FN
        WK(1) = 1.0D0 - XX*XX
        ABW2 = DABS(WK(1))
        WK(2) = DSQRT(ABW2)
        WK(7) = FN**CON2
        IF (ABW2.GT.0.27750D0) GO TO 80
C
C     ASYMPTOTIC EXPANSION
C     CASES NEAR X=FN, DABS(1.-(X/FN)**2).LE.0.2775
C     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES
C
C     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES
C
C     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)
C
        SA = 0.0D0
        IF (ABW2.EQ.0.0D0) GO TO 10
        SA = TOLS/DLOG(ABW2)
   10   SB = SA
        DO 20 I=1,5
          AKM = DMAX1(SA,2.0D0)
          KMAX(I) = INT(SNGL(AKM))
          SA = SA + SB
   20   CONTINUE
        KB = KMAX(5)
        KLAST = KB - 1
        SA = GAMA(KB)
        DO 30 K=1,KLAST
          KB = KB - 1
          SA = SA*WK(1) + GAMA(KB)
   30   CONTINUE
        Z = WK(1)*SA
        AZ = DABS(Z)
        RTZ = DSQRT(AZ)
        WK(3) = CON1*AZ*RTZ
        WK(4) = WK(3)*FN
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        IF(Z.LE.0.0D0) GO TO 35
        IF(WK(4).GT.ELIM) GO TO 75
        WK(6) = -WK(6)
   35   CONTINUE
        PHI = DSQRT(DSQRT(SA+SA+SA+SA))
C
C     B(ZETA) FOR S=0
C
        KB = KMAX(5)
        KLAST = KB - 1
        SB = BETA(KB,1)
        DO 40 K=1,KLAST
          KB = KB - 1
          SB = SB*WK(1) + BETA(KB,1)
   40   CONTINUE
        KSP1 = 1
        FN2 = FN*FN
        RFN2 = 1.0D0/FN2
        RDEN = 1.0D0
        ASUM = 1.0D0
        RELB = TOL*DABS(SB)
        BSUM = SB
        DO 60 KS=1,4
          KSP1 = KSP1 + 1
          RDEN = RDEN*RFN2
C
C     A(ZETA) AND B(ZETA) FOR S=1,2,3,4
C
          KSTEMP = 5 - KS
          KB = KMAX(KSTEMP)
          KLAST = KB - 1
          SA = ALFA(KB,KS)
          SB = BETA(KB,KSP1)
          DO 50 K=1,KLAST
            KB = KB - 1
            SA = SA*WK(1) + ALFA(KB,KS)
            SB = SB*WK(1) + BETA(KB,KSP1)
   50     CONTINUE
          TA = SA*RDEN
          TB = SB*RDEN
          ASUM = ASUM + TA
          BSUM = BSUM + TB
          IF (DABS(TA).LE.TOL .AND. DABS(TB).LE.RELB) GO TO 70
   60   CONTINUE
   70   CONTINUE
        BSUM = BSUM/(FN*WK(7))
        GO TO 160
C
   75   CONTINUE
        IFLW = 1
        RETURN
C
   80   CONTINUE
        UPOL(1) = 1.0D0
        TAU = 1.0D0/WK(2)
        T2 = 1.0D0/WK(1)
        IF (WK(1).GE.0.0D0) GO TO 90
C
C     CASES FOR (X/FN).GT.DSQRT(1.2775)
C
        WK(3) = DABS(WK(2)-DATAN(WK(2)))
        WK(4) = WK(3)*FN
        RCZ = -CON1/WK(4)
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        GO TO 100
   90   CONTINUE
C
C     CASES FOR (X/FN).LT.DSQRT(0.7225)
C
        WK(3) = DABS(DLOG((1.0D0+WK(2))/XX)-WK(2))
        WK(4) = WK(3)*FN
        RCZ = CON1/WK(4)
        IF(WK(4).GT.ELIM) GO TO 75
        Z32 = 1.5D0*WK(3)
        RTZ = Z32**CON2
        WK(7) = FN**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = WK(5)*WK(5)
  100   CONTINUE
        PHI = DSQRT((RTZ+RTZ)*TAU)
        TB = 1.0D0
        ASUM = 1.0D0
        TFN = TAU/FN
        RDEN=1.0D0/FN
        RFN2=RDEN*RDEN
        RDEN=1.0D0
        UPOL(2) = (C(1)*T2+C(2))*TFN
        CRZ32 = CON548*RCZ
        BSUM = UPOL(2) + CRZ32
        RELB = TOL*DABS(BSUM)
        AP = TFN
        KS = 0
        KP1 = 2
        RZDEN = RCZ
        L = 2
        ISETA=0
        ISETB=0
        DO 140 LR=2,8,2
C
C     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)
C
          LRP1 = LR + 1
          DO 120 K=LR,LRP1
            KS = KS + 1
            KP1 = KP1 + 1
            L = L + 1
            S1 = C(L)
            DO 110 J=2,KP1
              L = L + 1
              S1 = S1*T2 + C(L)
  110       CONTINUE
            AP = AP*TFN
            UPOL(KP1) = AP*S1
            CR(KS) = BR(KS)*RZDEN
            RZDEN = RZDEN*RCZ
            DR(KS) = AR(KS)*RZDEN
  120     CONTINUE
          SUMA = UPOL(LRP1)
          SUMB = UPOL(LR+2) + UPOL(LRP1)*CRZ32
          JU = LRP1
          DO 130 JR=1,LR
            JU = JU - 1
            SUMA = SUMA + CR(JR)*UPOL(JU)
            SUMB = SUMB + DR(JR)*UPOL(JU)
  130     CONTINUE
          RDEN=RDEN*RFN2
          TB = -TB
          IF (WK(1).GT.0.0D0) TB = DABS(TB)
          IF(RDEN.LT.TOL) GO TO 131
          ASUM = ASUM + SUMA*TB
          BSUM = BSUM + SUMB*TB
          GO TO 140
  131     IF(ISETA.EQ.1) GO TO 132
          IF(DABS(SUMA).LT.TOL) ISETA=1
          ASUM=ASUM+SUMA*TB
  132     IF(ISETB.EQ.1) GO TO 133
          IF(DABS(SUMB).LT.RELB) ISETB=1
          BSUM=BSUM+SUMB*TB
  133     IF(ISETA.EQ.1 .AND. ISETB.EQ.1) GO TO 150
  140   CONTINUE
  150   TB = WK(5)
        IF (WK(1).GT.0.0D0) TB = -TB
        BSUM = BSUM/TB
C
  160   CONTINUE
        CALL FUNJY(WK(6), WK(5), WK(4), FI, DFI)
        TA=1.0D0/TOL
        TB=D1MACH(1)*TA*1.0D+3
        IF(DABS(FI).GT.TB) GO TO 165
        FI=FI*TA
        DFI=DFI*TA
        PHI=PHI*TOL
  165   CONTINUE
        Y(JN) = FLGJY*PHI*(FI*ASUM+DFI*BSUM)/WK(7)
        FN = FN - FLGJY
  170 CONTINUE
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SAXPY-S DAXPY-D CAXPY-C),
C             LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DAXPY
C
      DOUBLE PRECISION DX(1),DY(1),DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBESI0(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DBESI0
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI0-S DBESI0-D),BESSEL FUNCTION,
C             FIRST KIND,HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION,ORDER ZERO,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. hyperbolic Bessel function of the first
C            kind of order zero.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C DBESI0(X) calculates the double precision modified (hyperbolic)
C Bessel function of the first kind of order zero and double
C precision argument X.
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DBSI0E,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DBESI0
      DOUBLE PRECISION X, BI0CS(18), XMAX, XSML, Y, D1MACH,
     1  DCSEVL, DBSI0E
      SAVE BI0 CS, NTI0, XSML, XMAX
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA NTI0, XSML, XMAX / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBESI0
      IF (NTI0.NE.0) GO TO 10
      NTI0 = INITDS (BI0CS, 18, 0.1*SNGL(D1MACH(3)))
      XSML = DSQRT (8.0D0*D1MACH(3))
      XMAX = DLOG (D1MACH(2))
C
 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20
C
      DBESI0 = 1.0D0
      IF (Y.GT.XSML) DBESI0 = 2.75D0 + DCSEVL (Y*Y/4.5D0-1.D0, BI0CS,
     1  NTI0)
      RETURN

 20   IF (Y.GT.XMAX) CALL XERROR ( 'DBESI0  DABS(X) SO BIG I0 OVERFLOWS'
     1, 35, 2, 2)

      DBESI0 = DEXP(Y) * DBSI0E(X)

      RETURN
      END
      SUBROUTINE DBESJ(X,ALPHA,N,Y,NZ)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DBESJ
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10A3
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(BESJ-S DBESJ-D),
C             BESSEL FUNCTION,J BESSEL FUNCTION,SPECIAL FUNCTIONS
C***AUTHOR  AMOS, D. E., (SNLA)
C           DANIEL, S. L., (SNLA)
C           WESTON, M. K., (SNLA)
C***PURPOSE  Compute an N member sequence of J Bessel functions
C            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA
C            and X.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C     Written by D. E. Amos, S. L. Daniel and M. K. Weston, January 1975
C
C     References
C         SAND-75-0147
C
C         CDC 6600 Subroutines IBESS and JBESS for Bessel Functions
C         I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0  by D. E. Amos, S. L.
C         Daniel, M. K. Weston. ACM Trans Math Software,3,pp 76-92
C         (1977)
C
C         Tables of Bessel Functions of Moderate or Large Orders,
C         NPL Mathematical Tables, Vol. 6, by F. W. J. Olver, Her
C         Majesty's Stationery Office, London, 1962.
C
C     Abstract  **** a double precision routine ****
C         DBESJ computes an N member sequence of J Bessel functions
C         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
C         A combination of the power series, the asymptotic expansion
C         for X to infinity and the uniform asymptotic expansion for
C         NU to infinity are applied over subdivisions of the (NU,X)
C         plane.  For values of (NU,X) not covered by one of these
C         formulae, the order is incremented or decremented by integer
C         values into a region where one of the formulae apply. Backward
C         recursion is applied to reduce orders by integer values except
C         where the entire sequence lies in the oscillatory region.  In
C         this case forward recursion is stable and values from the
C         asymptotic expansion for X to infinity start the recursion
C         when it is efficient to do so. Leading terms of the series and
C         uniform expansion are tested for underflow.  If a sequence is
C         requested and the last member would underflow, the result is
C         set to zero and the next lower order tried, etc., until a
C         member comes on scale or all members are set to zero.
C         Overflow cannot occur.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C         DBESJ calls DASYJY, DJAIRY, DLNGAM, D1MACH, I1MACH, XERROR
C
C     Description of Arguments
C
C         Input      X,ALPHA are double precision
C           X      - X .GE. 0.0D0
C           ALPHA  - order of first member of the sequence,
C                    ALPHA .GE. 0.0D0
C           N      - number of members in the sequence, N .GE. 1
C
C         Output     Y is double precision
C           Y      - a vector whose first N components contain
C                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
C           NZ     - number of components of Y set to zero due to
C                    underflow,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, last NZ components of Y set to zero,
C                             Y(K)=0.0D0, K=N-NZ+1,...,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Underflow  - a non-fatal error (NZ .NE. 0)
C***REFERENCES  CDC 6600 SUBROUTINES IBESS AND JBESS FOR BESSEL
C                 FUNCTIONS I(NU,X) AND J(NU,X), X .GE. 0, NU .GE. 0,
C                 BY D. E. AMOS, S. L.DANIEL, M. K. WESTON,  ACM
C                 TRANSACTIONS ON MATHEMATICALSOFTWARE, VOL. 3,
C                 PP. 76-92 (1977).
C***ROUTINES CALLED  D1MACH,DASYJY,DJAIRY,DLNGAM,I1MACH,XERROR
C***END PROLOGUE  DBESJ
      EXTERNAL DJAIRY
      INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN,
     1        NS,NZ
      INTEGER I1MACH
      DOUBLE PRECISION AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,
     1           EARG,ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,
     2           FNULIM,GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN,
     3           S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL,
     4           TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y,SLIM,RTOL
      SAVE RTWO, PDF, RTTP, PIDT, PP, INLIM, FNULIM
      DOUBLE PRECISION D1MACH, DLNGAM
      DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
      DATA RTWO,PDF,RTTP,PIDT                    / 1.34839972492648D+00,
     1 7.85398163397448D-01, 7.97884560802865D-01, 1.57079632679490D+00/
      DATA  PP(1),  PP(2),  PP(3),  PP(4)        / 8.72909153935547D+00,
     1 2.65693932265030D-01, 1.24578576865586D-01, 7.70133747430388D-04/
      DATA INLIM           /      150            /
      DATA FNULIM(1), FNULIM(2) /      100.0D0,     60.0D0     /
C***FIRST EXECUTABLE STATEMENT  DBESJ
      NZ = 0
      KT = 1
      NS=0
C     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
C     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
      TA = D1MACH(3)
      TOL = DMAX1(TA,1.0D-15)
      I1 = I1MACH(14) + 1
      I2 = I1MACH(15)
      TB = D1MACH(5)
      ELIM1 = 2.303D0*(DBLE(FLOAT(-I2))*TB-3.0D0)
      RTOL=1.0D0/TOL
      SLIM=D1MACH(1)*RTOL*1.0D+3
C     TOLLN = -LN(TOL)
      TOLLN = 2.303D0*TB*DBLE(FLOAT(I1))
      TOLLN = DMIN1(TOLLN,34.5388D0)
      IF (N-1) 720, 10, 20
   10 KT = 2
   20 NN = N
      IF (X) 730, 30, 80
   30 IF (ALPHA) 710, 40, 50
   40 Y(1) = 1.0D0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0D0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0D0) GO TO 710
C
      IALP = INT(SNGL(ALPHA))
      FNI = DBLE(FLOAT(IALP+N-1))
      FNF = ALPHA - DBLE(FLOAT(IALP))
      DFN = FNI + FNF
      FNU = DFN
      XO2 = X*0.5D0
      SXO2 = XO2*XO2
C
C     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
C     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
C     APPLIED.
C
      IF (SXO2.LE.(FNU+1.0D0)) GO TO 90
      TA = DMAX1(20.0D0,FNU)
      IF (X.GT.TA) GO TO 120
      IF (X.GT.12.0D0) GO TO 110
      XO2L = DLOG(XO2)
      NS = INT(SNGL(SXO2-FNU)) + 1
      GO TO 100
   90 FN = FNU
      FNP1 = FN + 1.0D0
      XO2L = DLOG(XO2)
      IS = KT
      IF (X.LE.0.50D0) GO TO 330
      NS = 0
  100 FNI = FNI + DBLE(FLOAT(NS))
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0D0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 330
  110 ANS = DMAX1(36.0D0-FNU,0.0D0)
      NS = INT(SNGL(ANS))
      FNI = FNI + DBLE(FLOAT(NS))
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 130
  120 CONTINUE
      RTX = DSQRT(X)
      TAU = RTWO*RTX
      TA = TAU + FNULIM(KT)
      IF (FNU.LE.TA) GO TO 480
      FN = FNU
      IS = KT
C
C     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
C
  130 CONTINUE
      I1 = IABS(3-IS)
      I1 = MAX0(I1,1)
      FLGJY = 1.0D0
      CALL DASYJY(DJAIRY,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
      IF(IFLW.NE.0) GO TO 380
      GO TO (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
      KT = 1
  320 IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 450
      GO TO 130
C
C     SERIES FOR (X/2)**2.LE.NU+1
C
  330 CONTINUE
      GLN = DLNGAM(FNP1)
      ARG = FN*XO2L - GLN
      IF (ARG.LT.(-ELIM1)) GO TO 400
      EARG = DEXP(ARG)
  340 CONTINUE
      S = 1.0D0
      IF (X.LT.TOL) GO TO 360
      AK = 3.0D0
      T2 = 1.0D0
      T = 1.0D0
      S1 = FN
      DO 350 K=1,17
        S2 = T2 + S1
        T = -T*SXO2/S2
        S = S + T
        IF (DABS(T).LT.TOL) GO TO 360
        T2 = T2 + AK
        AK = AK + 2.0D0
        S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 340
C
C     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
C     UNDERFLOW CAN ONLY OCCRFOR NS=0 SINCE THE ORDER MUST BE LARGER
C     THAN 36. THEREFORE, NS NEE NOT BE TESTED.
C
  380 Y(NN) = 0.0D0
      NN = NN - 1
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 390, 130
  390 KT = 2
      IS = 2
      GO TO 130
  400 Y(NN) = 0.0D0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 440, 410, 420
  410 KT = 2
      IS = 2
  420 IF (SXO2.LE.FNP1) GO TO 430
      GO TO 130
  430 ARG = ARG - XO2L + DLOG(FNP1)
      IF (ARG.LT.(-ELIM1)) GO TO 400
      GO TO 330
  440 NZ = N - NN
      RETURN
C
C     BACKWARD RECURSION SECTION
C
  450 CONTINUE
      IF(NS.NE.0) GO TO 451
      NZ = N - NN
      IF (KT.EQ.2) GO TO 470
C     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = TEMP(1)
      Y(NN-1) = TEMP(2)
      IF (NN.EQ.2) RETURN
  451 CONTINUE
      TRX = 2.0D0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      AK=1.0D0
      TA=TEMP(1)
      TB=TEMP(2)
      IF(DABS(TA).GT.SLIM) GO TO 455
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  455 CONTINUE
      KK=2
      IN=NS-1
      IF(IN.EQ.0) GO TO 690
      IF(NS.NE.0) GO TO 670
      K=NN-2
      DO 460 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  460 CONTINUE
      RETURN
  470 Y(1) = TEMP(2)
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
C     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
C     OF THE SEQUENCE IS ALSO IN THE REGION.
C
  480 CONTINUE
      IN = INT(SNGL(ALPHA-TAU+2.0D0))
      IF (IN.LE.0) GO TO 490
      IDALP = IALP - IN - 1
      KT = 1
      GO TO 500
  490 CONTINUE
      IDALP = IALP
      IN = 0
  500 IS = KT
      FIDAL = DBLE(FLOAT(IDALP))
      DALPHA = FIDAL + FNF
      ARG = X - PIDT*DALPHA - PDF
      SA = DSIN(ARG)
      SB = DCOS(ARG)
      COEF = RTTP/RTX
      ETX = 8.0D0*X
  510 CONTINUE
      DTM = FIDAL + FIDAL
      DTM = DTM*DTM
      TM = 0.0D0
      IF (FIDAL.EQ.0.0D0 .AND. DABS(FNF).LT.TOL) GO TO 520
      TM = 4.0D0*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
      TRX = DTM - 1.0D0
      T2 = (TRX+TM)/ETX
      S2 = T2
      RELB = TOL*DABS(T2)
      T1 = ETX
      S1 = 1.0D0
      FN = 1.0D0
      AK = 8.0D0
      DO 530 K=1,13
        T1 = T1 + ETX
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = -T2*AP/T1
        S1 = S1 + T2
        T1 = T1 + ETX
        AK = AK + 8.0D0
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = T2*AP/T1
        S2 = S2 + T2
        IF (DABS(T2).LE.RELB) GO TO 540
        AK = AK + 8.0D0
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
      IF(IS.EQ.2) GO TO 560
      FIDAL = FIDAL + 1.0D0
      DALPHA = FIDAL + FNF
      IS = 2
      TB = SA
      SA = -SB
      SB = TB
      GO TO 510
C
C     FORWARD RECURSION SECTION
C
  560 IF (KT.EQ.2) GO TO 470
      S1 = TEMP(1)
      S2 = TEMP(2)
      TX = 2.0D0/X
      TM = DALPHA*TX
      IF (IN.EQ.0) GO TO 580
C
C     FORWARD RECUR TO INDEX ALPHA
C
      DO 570 I=1,IN
        S = S2
        S2 = TM*S2 - S1
        TM = TM + TX
        S1 = S
  570 CONTINUE
      IF (NN.EQ.1) GO TO 600
      S = S2
      S2 = TM*S2 - S1
      TM = TM + TX
      S1 = S
  580 CONTINUE
C
C     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
C
      Y(1) = S1
      Y(2) = S2
      IF (NN.EQ.2) RETURN
      DO 590 I=3,NN
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TX
  590 CONTINUE
      RETURN
  600 Y(1) = S2
      RETURN
C
C     BACKWARD RECURSION WITH NORMALIZATION BY
C     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
C
  610 CONTINUE
C     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = DMAX1(3.0D0-FN,0.0D0)
      KM = INT(SNGL(AKM))
      TFN = FN + DBLE(FLOAT(KM))
      TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
      TA = XO2L - TA
      TB = -(1.0D0-1.5D0/TFN)/TFN
      AKM = TOLLN/(-TA+DSQRT(TA*TA-TOLLN*TB)) + 1.5D0
      IN = KM + INT(SNGL(AKM))
      GO TO 660
  620 CONTINUE
C     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      GLN = WK(3) + WK(2)
      IF (WK(6).GT.30.0D0) GO TO 640
      RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.0D0
      RZDEN = PP(1) + PP(2)*WK(6)
      TA = RZDEN/RDEN
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
      GO TO 650
  630 TB=(1.259921049D0+(0.1679894730D0+0.0887944358D0*WK(1))*WK(1))
     1 /WK(7)
      GO TO 650
  640 CONTINUE
      TA = 0.5D0*TOLLN/WK(4)
      TA=((0.0493827160D0*TA-0.1111111111D0)*TA+0.6666666667D0)*TA*WK(6)
      IF (WK(1).LT.0.10D0) GO TO 630
      TB = GLN/WK(5)
  650 IN = INT(SNGL(TA/TB+1.5D0))
      IF (IN.GT.INLIM) GO TO 310
  660 CONTINUE
      DTM = FNI + DBLE(FLOAT(IN))
      TRX = 2.0D0/X
      TM = (DTM+FNF)*TRX
      TA = 0.0D0
      TB = TOL
      KK = 1
      AK=1.0D0
  670 CONTINUE
C
C     BACKWARD RECUR UNINDEXED
C
      DO 680 I=1,IN
        S = TB
        TB = TM*TB - TA
        TA = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  680 CONTINUE
C     NORMALIZATION
      IF (KK.NE.1) GO TO 690
      S=TEMP(3)
      SA=TA/TB
      TA=S
      TB=S
      IF(DABS(S).GT.SLIM) GO TO 685
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  685 CONTINUE
      TA=TA*SA
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 670
  690 Y(NN) = TB*AK
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      K = NN - 1
      S=TB
      TB = TM*TB - TA
      TA=S
      Y(K)=TB*AK
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0D0
      TM = (DTM+FNF)*TRX
      K=NN-2
C
C     BACKWARD RECUR INDEXED
C
      DO 700 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  700 CONTINUE
      RETURN
C
C
C
  710 CONTINUE
      CALL XERROR('DBESJ - ORDER, ALPHA, LESS THAN ZERO.', 37, 2, 1)
      RETURN
  720 CONTINUE
      CALL XERROR('DBESJ - N LESS THAN ONE.', 24, 2, 1)
      RETURN
  730 CONTINUE
      CALL XERROR('DBESJ - X LESS THAN ZERO.', 25, 2, 1)
      RETURN
      END
      SUBROUTINE DBP (N,B,X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DBP
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C    Abstract
C       Computes values (at X) of the N+1 Bernstein basis functions 
C         of degree N on [0,1]
C
C    Description of Parameters
C       Input:  N (integer)           .GE. 0
C               X (double precision)  Abscissa for evaluation of polynomials
C       Output: B(0:N) (double precision array)
C                       B(I)= N!/(I!*(N-I)!) * (1.-X)**(N-I) * X**I
C                                       I=0,1,...,N
C***END PROLOGUE  DBP
C
      INTEGER  N,C,R
      DOUBLE PRECISION B(0:N),X

      IF (N .EQ. 0) THEN
        B(0) = 1.0D0
      ELSE IF (N .GT. 0) THEN
        DO 2 C=1,N
           IF (C .EQ. 1) THEN
                B(1) = X
           ELSE
                B(C) = X*B(C-1)
           END IF

           DO 1 R=C-1,1,-1
                B(R) = X*B(R-1) + (1.0D0-X)*B(R)
    1      CONTINUE

           IF (C .EQ. 1) THEN
                B(0) = 1.0D0-X
           ELSE
                B(0) = (1.0D0-X)*B(0)
           END IF
    2   CONTINUE
      END IF
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBSI0E(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DBSI0E
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C10B1
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(BESI0E-S DBSI0E-D),BESSEL FUNCTION,
C             EXPONENTIALLY SCALED,FIRST KIND,
C             HYPERBOLIC BESSEL FUNCTION,MODIFIED BESSEL FUNCTION,
C             ORDER ZERO,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. exponentially scaled hyperbolic Bessel
C            function of the first kind of order zero.
C***DESCRIPTION
C
C DBSI0E(X) calculates the double precision exponentially scaled
C modified (hyperbolic) Bessel function of the first kind of order
C zero for double precision argument X.  The result is the Bessel
C function I0(X) multiplied by EXP(-ABS(X)).
C
C Series for BI0        on the interval  0.          to  9.00000E+00
C                                        with weighted error   9.51E-34
C                                         log weighted error  33.02
C                               significant figures required  33.31
C                                    decimal places required  33.65
C
C Series for AI0        on the interval  1.25000E-01 to  3.33333E-01
C                                        with weighted error   2.74E-32
C                                         log weighted error  31.56
C                               significant figures required  30.15
C                                    decimal places required  32.39
C
C Series for AI02       on the interval  0.          to  1.25000E-01
C                                        with weighted error   1.97E-32
C                                         log weighted error  31.71
C                               significant figures required  30.15
C                                    decimal places required  32.63
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS
C***END PROLOGUE  DBSI0E
      DOUBLE PRECISION X, BI0CS(18), AI0CS(46), AI02CS(69),
     1  XSML, Y, D1MACH, DCSEVL
      SAVE BI0 CS, AI0 CS, AI02CS, NTI0, NTAI0, NTAI02, XSML
      DATA BI0 CS(  1) / -.7660547252 8391449510 8189497624 3285 D-1   /
      DATA BI0 CS(  2) / +.1927337953 9938082699 5240875088 1196 D+1   /
      DATA BI0 CS(  3) / +.2282644586 9203013389 3702929233 0415 D+0   /
      DATA BI0 CS(  4) / +.1304891466 7072904280 7933421069 1888 D-1   /
      DATA BI0 CS(  5) / +.4344270900 8164874513 7868268102 6107 D-3   /
      DATA BI0 CS(  6) / +.9422657686 0019346639 2317174411 8766 D-5   /
      DATA BI0 CS(  7) / +.1434006289 5106910799 6209187817 9957 D-6   /
      DATA BI0 CS(  8) / +.1613849069 6617490699 1541971999 4611 D-8   /
      DATA BI0 CS(  9) / +.1396650044 5356696994 9509270814 2522 D-10  /
      DATA BI0 CS( 10) / +.9579451725 5054453446 2752317189 3333 D-13  /
      DATA BI0 CS( 11) / +.5333981859 8625021310 1510774400 0000 D-15  /
      DATA BI0 CS( 12) / +.2458716088 4374707746 9678591999 9999 D-17  /
      DATA BI0 CS( 13) / +.9535680890 2487700269 4434133333 3333 D-20  /
      DATA BI0 CS( 14) / +.3154382039 7214273367 8933333333 3333 D-22  /
      DATA BI0 CS( 15) / +.9004564101 0946374314 6666666666 6666 D-25  /
      DATA BI0 CS( 16) / +.2240647369 1236700160 0000000000 0000 D-27  /
      DATA BI0 CS( 17) / +.4903034603 2428373333 3333333333 3333 D-30  /
      DATA BI0 CS( 18) / +.9508172606 1226666666 6666666666 6666 D-33  /
      DATA AI0 CS(  1) / +.7575994494 0237959427 2987203743 8 D-1      /
      DATA AI0 CS(  2) / +.7591380810 8233455072 9297873320 4 D-2      /
      DATA AI0 CS(  3) / +.4153131338 9237505018 6319749138 2 D-3      /
      DATA AI0 CS(  4) / +.1070076463 4390730735 8242970217 0 D-4      /
      DATA AI0 CS(  5) / -.7901179979 2128946607 5031948573 0 D-5      /
      DATA AI0 CS(  6) / -.7826143501 4387522697 8898980690 9 D-6      /
      DATA AI0 CS(  7) / +.2783849942 9488708063 8118538985 7 D-6      /
      DATA AI0 CS(  8) / +.8252472600 6120271919 6682913319 8 D-8      /
      DATA AI0 CS(  9) / -.1204463945 5201991790 5496089110 3 D-7      /
      DATA AI0 CS( 10) / +.1559648598 5060764436 1228752792 8 D-8      /
      DATA AI0 CS( 11) / +.2292556367 1033165434 7725480285 7 D-9      /
      DATA AI0 CS( 12) / -.1191622884 2790646036 7777423447 8 D-9      /
      DATA AI0 CS( 13) / +.1757854916 0324098302 1833124774 3 D-10     /
      DATA AI0 CS( 14) / +.1128224463 2189005171 4441135682 4 D-11     /
      DATA AI0 CS( 15) / -.1146848625 9272988777 2963387698 2 D-11     /
      DATA AI0 CS( 16) / +.2715592054 8036628726 4365192160 6 D-12     /
      DATA AI0 CS( 17) / -.2415874666 5626878384 4247572028 1 D-13     /
      DATA AI0 CS( 18) / -.6084469888 2551250646 0609963922 4 D-14     /
      DATA AI0 CS( 19) / +.3145705077 1754772937 0836026730 3 D-14     /
      DATA AI0 CS( 20) / -.7172212924 8711877179 6217505917 6 D-15     /
      DATA AI0 CS( 21) / +.7874493403 4541033960 8390960332 7 D-16     /
      DATA AI0 CS( 22) / +.1004802753 0094624023 4524457183 9 D-16     /
      DATA AI0 CS( 23) / -.7566895365 3505348534 2843588881 0 D-17     /
      DATA AI0 CS( 24) / +.2150380106 8761198878 1205128784 5 D-17     /
      DATA AI0 CS( 25) / -.3754858341 8308744291 5158445260 8 D-18     /
      DATA AI0 CS( 26) / +.2354065842 2269925769 0075710532 2 D-19     /
      DATA AI0 CS( 27) / +.1114667612 0479285302 2637335511 0 D-19     /
      DATA AI0 CS( 28) / -.5398891884 3969903786 9677932270 9 D-20     /
      DATA AI0 CS( 29) / +.1439598792 2407526770 4285840452 2 D-20     /
      DATA AI0 CS( 30) / -.2591916360 1110934064 6081840196 2 D-21     /
      DATA AI0 CS( 31) / +.2238133183 9985839074 3409229824 0 D-22     /
      DATA AI0 CS( 32) / +.5250672575 3647711727 7221683199 9 D-23     /
      DATA AI0 CS( 33) / -.3249904138 5332307841 7343228586 6 D-23     /
      DATA AI0 CS( 34) / +.9924214103 2050379278 5728471040 0 D-24     /
      DATA AI0 CS( 35) / -.2164992254 2446695231 4655429973 3 D-24     /
      DATA AI0 CS( 36) / +.3233609471 9435940839 7333299199 9 D-25     /
      DATA AI0 CS( 37) / -.1184620207 3967424898 2473386666 6 D-26     /
      DATA AI0 CS( 38) / -.1281671853 9504986505 4833868799 9 D-26     /
      DATA AI0 CS( 39) / +.5827015182 2793905116 0556885333 3 D-27     /
      DATA AI0 CS( 40) / -.1668222326 0261097193 6450150399 9 D-27     /
      DATA AI0 CS( 41) / +.3625309510 5415699757 0068480000 0 D-28     /
      DATA AI0 CS( 42) / -.5733627999 0557135899 4595839999 9 D-29     /
      DATA AI0 CS( 43) / +.3736796722 0630982296 4258133333 3 D-30     /
      DATA AI0 CS( 44) / +.1602073983 1568519633 6551253333 3 D-30     /
      DATA AI0 CS( 45) / -.8700424864 0572298845 2249599999 9 D-31     /
      DATA AI0 CS( 46) / +.2741320937 9374811456 0341333333 3 D-31     /
      DATA AI02CS(  1) / +.5449041101 4108831607 8960962268 0 D-1      /
      DATA AI02CS(  2) / +.3369116478 2556940898 9785662979 9 D-2      /
      DATA AI02CS(  3) / +.6889758346 9168239842 6263914301 1 D-4      /
      DATA AI02CS(  4) / +.2891370520 8347564829 6692402323 2 D-5      /
      DATA AI02CS(  5) / +.2048918589 4690637418 2760534093 1 D-6      /
      DATA AI02CS(  6) / +.2266668990 4981780645 9327743136 1 D-7      /
      DATA AI02CS(  7) / +.3396232025 7083863451 5084396952 3 D-8      /
      DATA AI02CS(  8) / +.4940602388 2249695891 0482449783 5 D-9      /
      DATA AI02CS(  9) / +.1188914710 7846438342 4084525196 3 D-10     /
      DATA AI02CS( 10) / -.3149916527 9632413645 3864862961 9 D-10     /
      DATA AI02CS( 11) / -.1321581184 0447713118 7540739926 7 D-10     /
      DATA AI02CS( 12) / -.1794178531 5068061177 7943574026 9 D-11     /
      DATA AI02CS( 13) / +.7180124451 3836662336 7106429346 9 D-12     /
      DATA AI02CS( 14) / +.3852778382 7421427011 4089801777 6 D-12     /
      DATA AI02CS( 15) / +.1540086217 5214098269 1325823339 7 D-13     /
      DATA AI02CS( 16) / -.4150569347 2872220866 2689972015 6 D-13     /
      DATA AI02CS( 17) / -.9554846698 8283076487 0214494312 5 D-14     /
      DATA AI02CS( 18) / +.3811680669 3526224207 4605535511 8 D-14     /
      DATA AI02CS( 19) / +.1772560133 0565263836 0493266675 8 D-14     /
      DATA AI02CS( 20) / -.3425485619 6772191346 1924790328 2 D-15     /
      DATA AI02CS( 21) / -.2827623980 5165834849 4205593759 4 D-15     /
      DATA AI02CS( 22) / +.3461222867 6974610930 9706250813 4 D-16     /
      DATA AI02CS( 23) / +.4465621420 2967599990 1042054284 3 D-16     /
      DATA AI02CS( 24) / -.4830504485 9441820712 5525403795 4 D-17     /
      DATA AI02CS( 25) / -.7233180487 8747539545 6227240924 5 D-17     /
      DATA AI02CS( 26) / +.9921475412 1736985988 8046093981 0 D-18     /
      DATA AI02CS( 27) / +.1193650890 8459820855 0439949924 2 D-17     /
      DATA AI02CS( 28) / -.2488709837 1508072357 2054491660 2 D-18     /
      DATA AI02CS( 29) / -.1938426454 1609059289 8469781132 6 D-18     /
      DATA AI02CS( 30) / +.6444656697 3734438687 8301949394 9 D-19     /
      DATA AI02CS( 31) / +.2886051596 2892243264 8171383073 4 D-19     /
      DATA AI02CS( 32) / -.1601954907 1749718070 6167156200 7 D-19     /
      DATA AI02CS( 33) / -.3270815010 5923147208 9193567485 9 D-20     /
      DATA AI02CS( 34) / +.3686932283 8264091811 4600723939 3 D-20     /
      DATA AI02CS( 35) / +.1268297648 0309501530 1359529710 9 D-22     /
      DATA AI02CS( 36) / -.7549825019 3772739076 9636664410 1 D-21     /
      DATA AI02CS( 37) / +.1502133571 3778353496 3712789053 4 D-21     /
      DATA AI02CS( 38) / +.1265195883 5096485349 3208799248 3 D-21     /
      DATA AI02CS( 39) / -.6100998370 0836807086 2940891600 2 D-22     /
      DATA AI02CS( 40) / -.1268809629 2601282643 6872095924 2 D-22     /
      DATA AI02CS( 41) / +.1661016099 8907414578 4038487490 5 D-22     /
      DATA AI02CS( 42) / -.1585194335 7658855793 7970504881 4 D-23     /
      DATA AI02CS( 43) / -.3302645405 9682178009 5381766755 6 D-23     /
      DATA AI02CS( 44) / +.1313580902 8392397817 4039623117 4 D-23     /
      DATA AI02CS( 45) / +.3689040246 6711567933 1425637280 4 D-24     /
      DATA AI02CS( 46) / -.4210141910 4616891492 1978247249 9 D-24     /
      DATA AI02CS( 47) / +.4791954591 0828657806 3171401373 0 D-25     /
      DATA AI02CS( 48) / +.8459470390 2218217952 9971707412 4 D-25     /
      DATA AI02CS( 49) / -.4039800940 8728324931 4607937181 0 D-25     /
      DATA AI02CS( 50) / -.6434714653 6504313473 0100850469 5 D-26     /
      DATA AI02CS( 51) / +.1225743398 8756659903 4464736990 5 D-25     /
      DATA AI02CS( 52) / -.2934391316 0257089231 9879821175 4 D-26     /
      DATA AI02CS( 53) / -.1961311309 1949829262 0371205728 9 D-26     /
      DATA AI02CS( 54) / +.1503520374 8221934241 6229900309 8 D-26     /
      DATA AI02CS( 55) / -.9588720515 7448265520 3386388206 9 D-28     /
      DATA AI02CS( 56) / -.3483339380 8170454863 9441108511 4 D-27     /
      DATA AI02CS( 57) / +.1690903610 2630436730 6244960725 6 D-27     /
      DATA AI02CS( 58) / +.1982866538 7356030438 9400115718 8 D-28     /
      DATA AI02CS( 59) / -.5317498081 4918162145 7583002528 4 D-28     /
      DATA AI02CS( 60) / +.1803306629 8883929462 3501450390 1 D-28     /
      DATA AI02CS( 61) / +.6213093341 4548931758 8405311242 2 D-29     /
      DATA AI02CS( 62) / -.7692189292 7721618632 0072806673 0 D-29     /
      DATA AI02CS( 63) / +.1858252826 1117025426 2556016596 3 D-29     /
      DATA AI02CS( 64) / +.1237585142 2813957248 9927154554 1 D-29     /
      DATA AI02CS( 65) / -.1102259120 4092238032 1779478779 2 D-29     /
      DATA AI02CS( 66) / +.1886287118 0397044900 7787447943 1 D-30     /
      DATA AI02CS( 67) / +.2160196872 2436589131 4903141406 0 D-30     /
      DATA AI02CS( 68) / -.1605454124 9197432005 8446594965 5 D-30     /
      DATA AI02CS( 69) / +.1965352984 5942906039 3884807331 8 D-31     /
      DATA NTI0, NTAI0, NTAI02, XSML / 3*0, 0.D0 /
C***FIRST EXECUTABLE STATEMENT  DBSI0E
      IF (NTI0.NE.0) GO TO 10
      ETA = 0.1*SNGL(D1MACH(3))
      NTI0 = INITDS (BI0CS, 18, ETA)
      NTAI0 = INITDS (AI0CS, 46, ETA)
      NTAI02 = INITDS (AI02CS, 69, ETA)
      XSML = DSQRT (8.0D0*D1MACH(3))

 10   Y = DABS(X)
      IF (Y.GT.3.0D0) GO TO 20

      DBSI0E = 1.0D0
      IF (Y.GT.XSML) DBSI0E = DEXP(-Y) * (2.75D0 +
     1  DCSEVL (Y*Y/4.5D0-1.D0, BI0CS, NTI0) )
      RETURN

 20   IF (Y.LE.8.D0) DBSI0E = (0.375D0 + DCSEVL ((48.D0/Y-11.D0)/5.D0,
     1  AI0CS, NTAI0))/DSQRT(Y)
      IF (Y.GT.8.D0) DBSI0E = (0.375D0 + DCSEVL (16.D0/Y-1.D0, AI02CS,
     1  NTAI02))/DSQRT(Y)

      RETURN
      END
      SUBROUTINE DCFFTB(N,C,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFFTB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A2
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Unnormalized inverse of DCFFTF.
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C  Subroutine DCFFTB computes the backward complex discrete Fourier
C  transform (the Fourier synthesis).  Equivalently, DCFFTB computes
C  a complex periodic sequence from its Fourier coefficients.
C  The transform is defined below at output parameter C.
C
C  A call of DCFFTF followed by a call of DCFFTB will multiply the
C  sequence by N.
C
C  The array WSAVE which is used by subroutine DCFFTB must be
C  initialized by calling subroutine DCFFTI(N,WSAVE).
C
C  Input Parameters
C
C
C  N      the length of the complex sequence C.  The method is
C         more efficient when N is the product of small primes.
C
C  C      a complex array of length N which contains the sequence
C
C  WSAVE   a d.p. work array which must be dimensioned at least 4*N+15
C          in the program that calls DCFFTB.  The WSAVE array must be
C          initialized by calling subroutine DCFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by DCFFTF and DCFFTB.
C
C  Output Parameters
C
C  C      For J=1,...,N
C
C             C(J)=the sum from K=1,...,N of
C
C                   C(K)*EXP(I*J*K*2*PI/N)
C
C                         where I=SQRT(-1)
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of subroutine DCFFTF or DCFFTB
C
C  *   References                                                      *
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCFTB1
C***END PROLOGUE  DCFFTB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       C(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DCFFTB
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE DCFFTF(N,C,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFFTF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A2
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Forward transform of a complex, periodic sequence.
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C  Subroutine DCFFTF computes the forward complex discrete Fourier
C  transform (the Fourier analysis).  Equivalently, DCFFTF computes
C  the Fourier coefficients of a complex periodic sequence.
C  The transform is defined below at output parameter C.
C
C  The transform is not normalized.  To obtain a normalized transform
C  the output must be divided by N.  Otherwise a call of DCFFTF
C  followed by a call of DCFFTB will multiply the sequence by N.
C
C  The array WSAVE which is used by subroutine DCFFTF must be
C  initialized by calling subroutine DCFFTI(N,WSAVE).
C
C  Input Parameters
C
C
C  N      the length of the complex sequence C.  The method is
C         more efficient when N is the product of small primes.
C
C  C      a complex array of length N which contains the sequence
C
C  WSAVE   a d.p. work array which must be dimensioned at least 4*N+15
C          in the program that calls DCFFTF.  The WSAVE array must be
C          initialized by calling subroutine DCFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by DCFFTF and DCFFTB.
C
C  Output Parameters
C
C  C      for J=1,...,N
C
C             C(J)=the sum from K=1,...,N of
C
C                   C(K)*EXP(-I*J*K*2*PI/N)
C
C                         where I=SQRT(-1)
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of subroutine DCFFTF or DCFFTB
C
C  *   References                                                      *
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCFTF1
C***END PROLOGUE  DCFFTF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       C(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DCFFTF
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE DCFFTI(N,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFFTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A2
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for DCFFTF and DCFFTB.
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C  Subroutine DCFFTI initializes the array WSAVE which is used in
C  both DCFFTF and DCFFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 4*N+15.
C          The same work array can be used for both DCFFTF and DCFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of DCFFTF or DCFFTB.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCFTI1
C***END PROLOGUE  DCFFTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DCFFTI
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
      SUBROUTINE DCFT2D(N,F,LDF,W,FORWD)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFT2D
C***DATE WRITTEN   870811   (YYMMDD)
C***REVISION DATE  870811   (YYMMDD)
C***CATEGORY NO.  J1B
C***KEYWORDS  TWO DIMENSIONAL FOURIER TRANSFORM, FFT
C***AUTHOR  KAHANER, DAVID K., (NBS)
C***PURPOSE  Two dimensional complex fast Fourier transform.
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C       Two dimensional fast Fourier transform, forward or backward
C          of complex N*N matrix F.
C
C   Input:
C     N: (INTEGER)        Number of rows and columns in the matrix F to be
C                            transformed. You must set N > 0, NOT CHECKED.
C     F: (COMPLEX*16)        Array of N*N complex values to be transformed.
C                            This array is overwritten on output.
C     LDF: (INTEGER)       Leading (first) dimension of the complex array F
C                             in the subroutine that calls DCFT2D. For example,
C                             if you declare F by either
C                               COMPLEX*16 F(0:15,0:20)  OR  F(16,21) then
C                             set LDF=16. You must have LDF >= N, NOT CHECKED.
C     W: (COMPLEX*16)         Array for internal use as work storage.
C                             Must be dimensioned by calling program to be
C                             at least 6N+15 D. P. words or 3N+8 COMPLEX*16 words.
C     FORWD: (LOGICAL)     Direction of transform. Set to .TRUE. for
C                            forward transform, set to .FALSE. for backward
C                            transform.
C
C
C   Output:
C     F: (COMPLEX*16)         Forward or reverse transformed input matrix.
C                          Output is unscaled, that is, a call to DCFT2D
C                          with FORWD=.TRUE. followed by a call to DCFT2D
C                          with FORWD=.FALSE. returns original data
C                          multiplied by N*N.
C
C
C   Remark:
C     For some applications it is desirable to have the transform scaled so
C        the center of the N by N frequency square corresponds to zero
C        frequency. The user can do this replacing the original input data
C        F(I,J) by F(I,J)*(-1.)**(I+J),  I,J =0,...,N-1.0D0
C
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED DCFFTI, DCFFTF, DCFFTB
C***END PROLOGUE  DCFT2D
      COMPLEX*16 F(0:LDF-1,0:*), W(0:*)
      LOGICAL FORWD
C
C  NOTE:  DOUBLE PRECISION COMPLEX IS NOT STANDARD FORTRAN.
C         SOME COMPILER VENDORS MAY FAIL TO IMPLEMENT IT
C         OR MAY IMPLEMENT IN A DIFFERENT MANNER.  IF YOUR
C         COMPILER DOES NOT ACCEPT THE COMPLEX*16 STATEMENT
C         WE RECOMMEND THAT YOU CHECK WITH YOUR VENDOR AND
C         MODIFY THE ABOVE ACCORDINGLY.
C
C***FIRST EXECUTABLE STATEMENT  DCFT2D
C   Find transform of each row, then of each column
C   First row transforms
      CALL DCFFTI(N,W(N))
      DO 10 I=0,N-1
C         Place row in beginning of array W
            DO 5 J=0,N-1
               W(J) = F(I,J)
    5       CONTINUE
C         Compute fft of row
            IF(FORWD) THEN
               CALL DCFFTF(N,W,W(N))
            ELSE
               CALL DCFFTB(N,W,W(N))
            ENDIF
C         Copy back to row of f
            DO 6 J=0,N-1
               F(I,J)=W(J)
    6       CONTINUE
   10 CONTINUE
C   Column transforms
      DO 20 J=0,N-1
C         Pass column of F to DCFFTF or DCFFTB by passing F(0,J).
C         Similarly DCFFTF or DCFFTB places results back there.
            IF(FORWD) THEN
               CALL DCFFTF(N, F(0,J), W(N))
            ELSE
               CALL DCFFTB(N, F(0,J), W(N))
            ENDIF
   20 CONTINUE
      RETURN
      END
      SUBROUTINE DCFTB1(N,C,CH,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFTB1
C***REFER TO  DCFFTB
C***ROUTINES CALLED  DPASSB,DPASB2,DPASB3,DPASB4,DPASB5
C***END PROLOGUE  DCFTB1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
C***FIRST EXECUTABLE STATEMENT  DCFTB1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE DCFTF1(N,C,CH,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFTF1
C***REFER TO  DCFFTF
C***ROUTINES CALLED  DPASSF,DPASF2,DPASF3,DPASF4,DPASF5
C***END PROLOGUE  DCFTF1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
C***FIRST EXECUTABLE STATEMENT  DCFTF1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASF2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE DCFTI1(N,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCFTI1
C***REFER TO  DCFFTI
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DCFTI1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
C***FIRST EXECUTABLE STATEMENT  DCFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0D0*ATAN(1.0D0)
      ARGH = TPI/DBLE(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.0D0
            WA(I) = 0.0D0
            LD = LD+L1
            FI = 0.0D0
            ARGLD = DBLE(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.0D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
      SUBROUTINE DCHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCHFDV
C***DATE WRITTEN   811019   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(CHFDV-S DCHFDV-D),
C             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
C             CUBIC POLYNOMIAL EVALUATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
C            first derivative at an array of points.  While designed for
C            use by DPCHFD, it may be useful directly as an evaluator
C            for a piecewise cubic Hermite function in applications,
C            such as graphing, where the interval is known in advance.
C            If only function values are required, use DCHFEV instead.
C***DESCRIPTION
C
C       **** Double Precision version of CHFDV ****
C
C        DCHFDV:  Cubic Hermite Function and Derivative Evaluator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
C     its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use DCHFEV, instead.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  NE, NEXT(2), IERR
C        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE),
C                          DE(NE)
C
C        CALL  DCHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
C
C   Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.EQ.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real*8 array of points at which the functions are to
C           be evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned in NEXT.
C
C     FE -- (output) real*8 array of values of the cubic function
C           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     DE -- (output) real*8 array of values of the first derivative of
C           the same function at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if NE.LT.1 .
C              IERR = -2  if X1.EQ.X2 .
C                (Output arrays have not been changed in either case.)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DCHFDV
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-03   Minor cosmetic changes for release 1.
C     87-07-07   Corrected XERROR calls for d.p. names(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFDV to CHFDV wherever it occurs,
C        b. Change the double precision declaration to real,
C        c. Change the constant Zero to single precision, and
C        d. Change the names of the Fortran functions:  AMAX1, AMIN1.
C
C  DECLARE ARGUMENTS.
C
      INTEGER NE, NEXT(2), IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I
      DOUBLE PRECISION  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X,
     *  XMI, XMA, ZERO
      DATA  ZERO /0.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DCHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = DMIN1(ZERO, H)
      XMA = DMAX1(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
      CALL XERROR ('DCHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
     *           , 51, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERROR ('DCHFDV -- INTERVAL ENDPOINTS EQUAL'
     *           , 34, IERR, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCHFIV(X1,X2,F1,F2,D1,D2,A,B,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCHFIV
C***REFER TO  DPCHIA
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DCHFIV:  Cubic Hermite Function Integral Evaluator.
C
C     Called by  DPCHIA  to evaluate the integral of a single cubic (in
C     Hermite form) over an arbitrary interval (A,B).
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  IERR
C        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
C        DOUBLE PRECISION  VALUE, DCHFIV
C
C        VALUE = DCHFIV (X1, X2, F1, F2, D1, D2, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     X1,X2 -- (input) endpoints if interval of definition of cubic.
C           (Must be distinct.  Error return if not.)
C
C     F1,F2 -- (input) function values at the ends of the interval.
C
C     D1,D2 -- (input) derivative values at the ends of the interval.
C
C     A,B -- (input) endpoints of interval of integration.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0 (no errors).
C           "Recoverable errors":
C              IERR = -1  if X1.EQ.X2 .
C                (VALUE has not been set in this case.)
C
C***END PROLOGUE  DCHFIV
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFIV to CHFIV wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants HALF, TWO, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
C
C  DECLARE LOCAL VARIABLES.
C
      DOUBLE PRECISION  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2,
     *      PHIB1, PHIB2, PSIA1, PSIA2, PSIB1, PSIB2, TA1, TA2, TB1,
     *      TB2, THREE, TWO, UA1, UA2, UB1, UB2
C
C  INITIALIZE.
C
      DATA  HALF/.5D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/, SIX/6.D0/
C
C  VALIDITY CHECK INPUT.
C
C***FIRST EXECUTABLE STATEMENT  DCHFIV
      IF (X1 .EQ. X2)  GO TO 5001
      IERR = 0
C
C  COMPUTE INTEGRAL.
C
      H = X2 - X1
      TA1 = (A - X1) / H
      TA2 = (X2 - A) / H
      TB1 = (B - X1) / H
      TB2 = (X2 - B) / H
C
      UA1 = TA1**3
      PHIA1 = UA1 * (TWO - TA1)
      PSIA1 = UA1 * (THREE*TA1 - FOUR)
      UA2 = TA2**3
      PHIA2 =  UA2 * (TWO - TA2)
      PSIA2 = -UA2 * (THREE*TA2 - FOUR)
C
      UB1 = TB1**3
      PHIB1 = UB1 * (TWO - TB1)
      PSIB1 = UB1 * (THREE*TB1 - FOUR)
      UB2 = TB2**3
      PHIB2 =  UB2 * (TWO - TB2)
      PSIB2 = -UB2 * (THREE*TB2 - FOUR)
C
      FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
      DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
C
C  RETURN VALUE.
C
      DCHFIV = (HALF*H) * (FTERM + DTERM)
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('DCHFIV -- X1 EQUAL TO X2'
     *           , 0, IERR, 1)
      RETURN
      END
      SUBROUTINE DCKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCKDER
C***DATE WRITTEN   800301   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  F3,G4C
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(CHKDER-S DCKDER-D),
C             GRADIENTS,JACOBIAN,MINPACK,NONLINEAR
C***AUTHOR  HIEBERT K.L. (SNLA)
C***PURPOSE  Checks the gradients of M nonlinear functions in N
C            variables, evaluated at a point X, for consistency
C            with the functions themselves.
C***DESCRIPTION
C
C   This subroutine is a companion routine to DNSQ and DNSQE. It may
C   be used to check the coding of the Jacobian calculation.
C
C     SUBROUTINE DCKDER
C
C     This subroutine checks the gradients of M nonlinear functions
C     in N variables, evaluated at a point X, for consistency with
C     the functions themselves. The user must call DCKDER twice,
C     first with MODE = 1 and then with MODE = 2.
C
C     MODE = 1. On input, X must contain the point of evaluation.
C               On output, XP is set to a neighboring point.
C
C     MODE = 2. On input, FVEC must contain the functions and the
C                         rows of FJAC must contain the gradients
C                         of the respective functions each evaluated
C                         at X, and FVECP must contain the functions
C                         evaluated at XP.
C               On output, ERR contains measures of correctness of
C                          the respective gradients.
C
C     The subroutine does not perform reliably if cancellation or
C     rounding errors cause a severe loss of significance in the
C     evaluation of a function. Therefore, none of the components
C     of X should be unusually small (in particular, zero) or any
C     other value which may cause loss of significance.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE DCKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of functions.
C
C       N is a positive integer input variable set to the number
C         of variables.
C
C       X is an input array of length N.
C
C       FVEC is an array of length M. On input when MODE = 2,
C         FVEC must contain the functions evaluated at X.
C
C       FJAC is an M by N array. On input when MODE = 2,
C         the rows of FJAC must contain the gradients of
C         the respective functions evaluated at X.
C
C       LDFJAC is a positive integer input parameter not less than M
C         which specifies the leading dimension of the array FJAC.
C
C       XP is an array of length N. On output when MODE = 1,
C         XP is set to a neighboring point of X.
C
C       FVECP is an array of length M. On input when MODE = 2,
C         FVECP must contain the functions evaluated at XP.
C
C       MODE is an integer input variable set to 1 on the first call
C         and 2 on the second. Other values of MODE are equivalent
C         to MODE = 1.
C
C       ERR is an array of length M. On output when MODE = 2,
C         ERR contains measures of correctness of the respective
C         gradients. If there is no severe loss of significance,
C         then if ERR(I) is 1.0 the I-th gradient is correct,
C         while if ERR(I) is 0.0 the I-th gradient is incorrect.
C         For values of ERR between 0.0 and 1.0, the categorization
C         is less certain. In general, a value of ERR(I) greater
C         than 0.5 indicates that the I-th gradient is probably
C         correct, while a value of ERR(I) less than 0.5 indicates
C         that the I-th gradient is probably incorrect.
C
C     Subprograms called
C
C       SLATEC supplied ... D1MACH
C
C       FORTRAN supplied ... DABS,DLOG10,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***REFERENCES  POWELL, M. J. D.
C               A HYBRID METHOD FOR NONLINEAR EQUATIONS.
C               NUMERICAL METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS,
C               P. RABINOWITZ, EDITOR.  GORDON AND BREACH, 1970.
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DCKDER
      INTEGER I, J, LDFJAC, M, MODE, N
      DOUBLE PRECISION D1MACH, EPS, EPSF, EPSLOG, EPSMCH, ERR(M),
     1     FACTOR, FJAC(LDFJAC,N), FVEC(M), FVECP(M), ONE, TEMP, X(N),
     2     XP(N), ZERO
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
C***FIRST EXECUTABLE STATEMENT  DCKDER
      EPSMCH = D1MACH(4)
C
      EPS = DSQRT(EPSMCH)
C
      IF (MODE .EQ. 2) GO TO 20
C
C        MODE = 1.
C
         DO 10 J = 1, N
            TEMP = EPS*DABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = EPS
            XP(J) = X(J) + TEMP
   10       CONTINUE
         GO TO 70
   20 CONTINUE
C
C        MODE = 2.
C
         EPSF = FACTOR*EPSMCH
         EPSLOG = DLOG10(EPS)
         DO 30 I = 1, M
            ERR(I) = ZERO
   30       CONTINUE
         DO 50 J = 1, N
            TEMP = DABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = ONE
            DO 40 I = 1, M
               ERR(I) = ERR(I) + TEMP*FJAC(I,J)
   40          CONTINUE
   50       CONTINUE
         DO 60 I = 1, M
            TEMP = ONE
            IF (FVEC(I) .NE. ZERO .AND. FVECP(I) .NE. ZERO
     1          .AND. DABS(FVECP(I)-FVEC(I)) .GE. EPSF*DABS(FVEC(I)))
     2         TEMP = EPS*DABS((FVECP(I)-FVEC(I))/EPS-ERR(I))
     3                /(DABS(FVEC(I)) + DABS(FVECP(I)))
            ERR(I) = ONE
            IF (TEMP .GT. EPSMCH .AND. TEMP .LT. EPS)
     1         ERR(I) = (DLOG10(TEMP) - EPSLOG)/EPSLOG
            IF (TEMP .GE. EPS) ERR(I) = ZERO
   60       CONTINUE
   70 CONTINUE
C
      RETURN
      END
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCOPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SCOPY-S DCOPY-D CCOPY-C),COPY,
C             LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector copy y = x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DCOPY
C
      DOUBLE PRECISION DX(1),DY(1)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS=N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DX(I)
   70     CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DCSEVL(X,A,N)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DCSEVL
C***DATE WRITTEN   770401   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,FNLIB,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Evaluate the double precision N-term Chebyshev series A
C            at X.
C***DESCRIPTION
C
C Evaluate the N-term Chebyshev series A at X.  Adapted from
C R. Broucke, Algorithm 446, C.A.C.M., 16, 254 (1973).
C W. Fullerton, C-3, Los Alamos Scientific Laboratory.
C
C       Input Arguments --
C X    double precision value at which the series is to be evaluated.
C A    double precision array of N terms of a Chebyshev series.  In
C      evaluating A, only half of the first coefficient is summed.
C N    number of terms in array A.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DCSEVL
C
       DOUBLE PRECISION A(N),X,TWOX,B0,B1,B2
C***FIRST EXECUTABLE STATEMENT  DCSEVL
       IF(N.LT.1)CALL XERROR( 'DCSEVL  NUMBER OF TERMS LE 0', 28, 2,2)
       IF(N.GT.1000) CALL XERROR ( 'DCSEVL  NUMBER OF TERMS GT 1000',
     1   31, 3, 2)
       IF ((X.LT.-1.D0) .OR. (X.GT.1.D0)) CALL XERROR ( 'DCSEVL  X OUTSI
     1DE (-1,+1)', 25, 1, 1)
C
       TWOX = 2.0D0*X
       B1 = 0.D0
       B0=0.D0
       DO 10 I=1,N
         B2=B1
         B1=B0
         NI = N - I + 1
         B0 = TWOX*B1 - B2 + A(NI)
 10    CONTINUE
C
       DCSEVL = 0.5D0 * (B0-B2)
C
       RETURN
      END
      SUBROUTINE DDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
     8   NDE,NQ,T,USERS,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D,JSTATE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDCOR
C***REFER TO  DDRIV3
C  Subroutine DDCOR is called to compute corrections to the Y array.
C  In the case of functional iteration, update Y directly from the
C  result of the last call to F.
C  In the case of the chord method, compute the corrector error and
C  solve the linear system with that as right hand side and DFDY as
C  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
C  or 5.
C***ROUTINES CALLED  DGESL,DGBSL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870401   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDCOR
      DOUBLE PRECISION A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,
     8     SAVE1(*), SAVE2(*), DNRM2, T, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL EVALFA
C***FIRST EXECUTABLE STATEMENT  DDCOR
      IF (MITER .EQ. 0) THEN
        DO 100 I = 1,N
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        D = DNRM2(N, SAVE1, 1)/SQRT(DBLE(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 150 I = 1,N
 150        SAVE2(I) = H*SAVE2(I)
          DO 160 J = 1,N
            DO 160 I = 1,N
 160          SAVE2(I) = SAVE2(I) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL DGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        DO 200 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 200      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 250 I = 1,N
 250        SAVE2(I) = H*SAVE2(I)
          MW = ML + 1 + MU
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
              I3 = I + J - MW
 260          SAVE2(I3) = SAVE2(I3) - A(I,J)*(YH(J,2) + SAVE1(J))
        ELSE IF (IMPL .EQ. 2) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL DGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        DO 300 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 300      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      END IF
      END
      SUBROUTINE DDCST (MAXORD,MINT,ISWFLG,EL,TQ)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDCST
C***REFER TO  DDRIV3
C  DDCST is called by DDNTL and sets coefficients used by the core
C  integrator DDSTP.  The array EL determines the basic method.
C  The array TQ is involved in adjusting the step size in relation
C  to truncation error.  EL and TQ depend upon MINT, and are calculated
C  for orders 1 to MAXORD(.LE. 12).  For each order NQ, the coefficients
C  EL are calculated from the generating polynomial:
C    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.
C  For the implicit Adams methods, L(T) is given by
C    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0,
C    where      K = factorial(NQ-1).
C  For the Gear methods,
C    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K,
C    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).
C  For each order NQ, there are three components of TQ.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870216   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDCST
      DOUBLE PRECISION EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
C***FIRST EXECUTABLE STATEMENT  DDCST
      FACTRL(1) = 1.D0
      DO 10 I = 2,MAXORD
 10     FACTRL(I) = DBLE(I)*FACTRL(I-1)
C                                             COMPUTE ADAMS COEFFICIENTS
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.D0
        DO 40 I = 1,MAXORD+1
          SUM = 0.D0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        EL(2,2) = 1.D0
        EL(3,2) = 1.D0
        DO 60 J = 3,MAXORD
          EL(2,J) = FACTRL(J-1)
          DO 50 I = 3,J
 50         EL(I,J) = DBLE(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60       EL(J+1,J) = 1.D0
        DO 80 J = 2,MAXORD
          EL(1,J) = EL(1,J-1) + GAMMA(J)
          EL(2,J) = 1.D0
          DO 80 I = 3,J+1
 80         EL(I,J) = EL(I,J)/(DBLE(I-1)*FACTRL(J-1))
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.D0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.D0/GAMMA(J+1)
 100      TQ(3,J) = -1.D0/GAMMA(J+2)
C                                              COMPUTE GEAR COEFFICIENTS
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        DO 130 J = 2,MAXORD
          EL(1,J) = FACTRL(J)
          DO 120 I = 2,J
 120        EL(I,J) = DBLE(J)*EL(I,J-1) + EL(I-1,J-1)
 130      EL(J+1,J) = 1.D0
        SUM = 1.D0
        DO 150 J = 2,MAXORD
          SUM = SUM + 1.D0/DBLE(J)
          DO 150 I = 1,J+1
 150        EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.D0/FACTRL(J-1)
          TQ(2,J) = DBLE(J+1)/EL(1,J)
 170      TQ(3,J) = DBLE(J+2)/EL(1,J)
      END IF
C                          Compute constants used in the stiffness test.
C                          These are the ratio of TQ(2,NQ) for the Gear
C                          methods to those for the Adams methods.
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.D0
          DO 190 I = 1,MXRD
            SUM = 0.D0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/DBLE(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        SUM = 1.D0
        DO 200 I = 2,MXRD
          SUM = SUM + 1.D0/DBLE(I)
 200      EL(1+I,1) = -DBLE(I+1)*SUM*GAMMA(I+1)
      END IF
      END
      SUBROUTINE DDIRECT (N,DATA,AZERO,A,B)

c*********************************************************************72
c
C Direct use of definitions to compute real DFT
C No simplifications...SLOW
c
      DOUBLE PRECISION DATA(0:*),A(*),B(*),AZERO,TPN
      AZERO = 0.0D0
      DO 1 J = 0,N-1
         AZERO = AZERO+DATA(J)
    1 CONTINUE
      AZERO = AZERO/N
      TPN = 2.0D0*ASIN(1.0D0)*2.0D0/N
      DO 20 K = 1,N/2
         A(K) = 0D0
         B(K) = 0D0
         DO 10 J = 0,N-1
            A(K) = A(K) + DATA(J)*COS(J*K*TPN)
            B(K) = B(K) + DATA(J)*SIN(J*K*TPN)
   10    CONTINUE
         IF (K .NE. N/2)THEN
            A(K) = A(K)*(2.0D0/N)
            B(K) = B(K)*(2.0D0/N)
         ELSE
            A(K) = A(K)*(1.0D0/N)
            B(K) = B(K)*(1.0D0/N)
         END IF
   20 CONTINUE
      RETURN
      END 
      SUBROUTINE DDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,UROUND,USERS,Y,YWT,H,MNTOLD,
     8   MTROLD,NFE,RC,YH,A,CONVRG,EL,FAC,IER,IPVT,NQ,NWAIT,RH,RMAX,
     8   SAVE2,TQ,TREND,ISWFLG,JSTATE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDNTL
C***REFER TO  DDRIV3
C  Subroutine DDNTL is called to set parameters on the first call
C  to DDSTP, on an internal restart, or when the user has altered
C  MINT, MITER, and/or H.
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, DDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C***ROUTINES CALLED  DDCST,DDSCL,DGEFA,DGESL,DGBFA,DGBSL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870810   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDNTL
      DOUBLE PRECISION A(MATDIM,*), EL(13,12), EPS, FAC(*), H, HMAX,
     8     HOLD, OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX,
     8     SMIN, DNRM2, SUM, SUM0, T, TQ(3,12), TREND, UROUND, Y(*),
     8     YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)
C***FIRST EXECUTABLE STATEMENT  DDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) THEN
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RMAX = RMINIT
        END IF
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          RETURN
        END IF
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
            IF (N .EQ. 0) THEN
              JSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              IF (N .EQ. 0) THEN
                JSTATE = 9
                RETURN
              END IF
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            IF (N .EQ. 0) THEN
              JSTATE = 9
              RETURN
            END IF
            DO 150 I = 1,NDE
              IF (A(I,1) .EQ. 0.D0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.D0
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = DNRM2(NDE, SAVE1, 1)
        SUM0 = 1.D0/MAX(1.D0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.D0 + (SMIN/SMAX)**2)/SQRT(DBLE(NDE))
        H = SIGN(MIN(2.D0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
        IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. ISWFLG .EQ. 3) THEN
          DO 20 I = 1,N
 20         FAC(I) = SQRT(UROUND)
        END IF
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.D0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          RH = H/HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  HOLD, RC, RH, YH)
        END IF
      END IF
      END
      SUBROUTINE DDNTP (H,K,N,NQ,T,TOUT,YH,Y)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDNTP
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C***REFER TO  DDRIV3
C   Subroutine DDNTP interpolates the K-th derivative of Y at TOUT,
C   using the data in the YH array.  If K has a value greater than NQ,
C   the NQ-th derivative is calculated.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870216   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDNTP
      DOUBLE PRECISION FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDNTP
      IF (K .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        KUSED = MIN(K, NQ)
        FACTOR = 1.D0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*DBLE(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        DO 80 JJ = KUSED+1,NQ
          J = K + 1 + NQ - JJ
          FACTOR = 1.D0
          DO 60 KK = 1,KUSED
 60         FACTOR = FACTOR*DBLE(J-KK)
          DO 70 I = 1,N
 70         Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80       CONTINUE
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
      SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDOGLG
C***REFER TO  DNSQ,DNSQE
C
C       *********
C
C     SUBROUTINE DDOGLG
C
C     Given an M by N matrix A, an N by N nonsingular diagonal
C     matrix D, an M-vector B, and a positive number DELTA, the
C     problem is to determine the convex combination X of the
C     Gauss-Newton and scaled gradient directions that minimizes
C     (A*X - B) in the least squares sense, subject to the
C     restriction that the Euclidean norm of D*X be at most DELTA.
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization of A. That is, if A = Q*R, where Q has
C     orthogonal columns and R is an upper triangular matrix,
C     then DDOGLG expects the full upper triangle of R and
C     the first N components of (Q transpose)*B.
C
C     The subroutine statement is
C
C       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     where
C
C       N is a positive integer input variable set to the order of R.
C
C       R is an input array of length LR which must contain the upper
C         triangular matrix R stored by rows.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       DIAG is an input array of length N which must contain the
C         diagonal elements of the matrix D.
C
C       QTB is an input array of length N which must contain the first
C         N elements of the vector (Q transpose)*B.
C
C       DELTA is a positive input variable which specifies an upper
C         bound on the Euclidean norm of D*X.
C
C       X is an output array of length N which contains the desired
C         convex combination of the Gauss-Newton direction and the
C         scaled gradient direction.
C
C       WA1 and WA2 are work arrays of length N.
C
C     Subprograms called
C
C       SLATEC-supplied ... D1MACH,DENORM
C
C       FORTRAN-supplied ... DABS,DMAX1,DMIN1,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***ROUTINES CALLED  D1MACH,DENORM
C***END PROLOGUE  DDOGLG
      DOUBLE PRECISION D1MACH,DENORM
      DOUBLE PRECISION DABS, DMAX1, DMIN1, DSQRT
      INTEGER I, J, JJ, JP1, K, L, LR, N
      DOUBLE PRECISION ALPHA, BNORM, DELTA, DIAG(N), EPSMCH, GNORM,
     1     ONE, QNORM, QTB(N), R(LR), SGNORM, SUM, TEMP, WA1(N),
     2     WA2(N), X(N), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
C***FIRST EXECUTABLE STATEMENT  DDOGLG
      EPSMCH = D1MACH(4)
C
C     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
C
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N .LT. JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (TEMP .NE. ZERO) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = DMAX1(TEMP,DABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
C
C     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
C
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
      QNORM = DENORM(N,WA2)
      IF (QNORM .LE. DELTA) GO TO 140
C
C     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
C     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
C
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
C
C     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
C     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
C
      GNORM = DENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM .EQ. ZERO) GO TO 120
C
C     CALCULATE THE POINT ALONG THE SCALED GRADIENT
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = DENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
C
C     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
C
      ALPHA = ZERO
      IF (SGNORM .GE. DELTA) GO TO 120
C
C     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
C     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      BNORM = DENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2
     1       + DSQRT((TEMP-(DELTA/QNORM))**2
     2               +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
C
C     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
C     DIRECTION AND THE SCALED GRADIENT DIRECTION.
C
      TEMP = (ONE - ALPHA)*DMIN1(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SDOT-S DDOT-D CDOTU-C),
C             INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. inner product of d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DDOT
C
      DOUBLE PRECISION DX(1),DY(1)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE DDPSC (KSGN,N,NQ,YH)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDPSC
C***REFER TO  DDRIV3
C     This subroutine computes the predicted YH values by effectively
C     multiplying the YH array by the Pascal triangle matrix when KSGN
C     is +1, and performs the inverse function when KSGN is -1.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDPSC
      DOUBLE PRECISION YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDPSC
      IF (KSGN .GT. 0) THEN
        DO 10 J1 = 1,NQ
          DO 10 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 10 I = 1,N
 10           YH(I,J) = YH(I,J) + YH(I,J+1)
      ELSE
        DO 30 J1 = 1,NQ
          DO 30 J2 = J1,NQ
            J = NQ - J2 + J1
            DO 30 I = 1,N
 30           YH(I,J) = YH(I,J) - YH(I,J+1)
      END IF
      END
      SUBROUTINE DDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
     8   NQ,SAVE2,T,USERS,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,FAC,IER,IPVT,
     8   SAVE1,ISWFLG,BND,JSTATE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDPST
C***REFER TO  DDRIV3
C  Subroutine DDPST is called to reevaluate the partials.
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C***ROUTINES CALLED  DGEFA,DGBFA,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870401   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDPST
      DOUBLE PRECISION A(MATDIM,*), BL, BND, BP, BR, BU, DFDY(MATDIM,*),
     8     DFDYMX, DIFF, DY, EL(13,12), FAC(*), FACMAX, FACMIN, FACTOR,
     8     H, SAVE1(*), SAVE2(*), SCALE, DNRM2, T, UROUND, Y(*),
     8     YH(N,*), YJ, YS, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER
      PARAMETER(FACMAX = .5D0)
C***FIRST EXECUTABLE STATEMENT  DDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BU = UROUND**(.25D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          DO 170 J = 1,N
            YS = MAX(ABS(YWT(J)), ABS(Y(J)))
 120        DY = FAC(J)*YS
            IF (DY .EQ. 0.D0) THEN
              IF (FAC(J) .LT. FACMAX) THEN
                FAC(J) = MIN(100.D0*FAC(J), FACMAX)
                GO TO 120
              ELSE
                DY = YS
              END IF
            END IF
            IF (NQ .EQ. 1) THEN
              DY = SIGN(DY, SAVE2(J))
            ELSE
              DY = SIGN(DY, YH(J,3))
            END IF
            DY = (Y(J) + DY) - Y(J)
            YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            Y(J) = YJ
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
C                                                                 Step 1
            DIFF = ABS(SAVE2(1) - SAVE1(1))
            IMAX = 1
            DO 150 I = 2,N
              IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                IMAX = I
                DIFF = ABS(SAVE2(I) - SAVE1(I))
              END IF
 150          CONTINUE
C                                                                 Step 2
            IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT. 0.D0) THEN
              SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
              IF (DIFF .GT. BU*SCALE) THEN
                FAC(J) = MAX(FACMIN, FAC(J)*.1D0)
              ELSE IF (BR*SCALE .LE. DIFF .AND. DIFF .LE. BL*SCALE) THEN
                FAC(J) = MIN(FAC(J)*10.D0, FACMAX)
C                                                                 Step 4
              ELSE IF (DIFF .LT. BR*SCALE) THEN
                FAC(J) = MIN(BP*FAC(J), FACMAX)
              END IF
            END IF
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
          NFE = NFE + N
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL DGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          IF (N .EQ. 0) THEN
            JSTATE = 8
            RETURN
          END IF
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          BR = UROUND**(.875D0)
          BL = UROUND**(.75D0)
          BU = UROUND**(.25D0)
          BP = UROUND**(-.15D0)
          FACMIN = UROUND**(.78D0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 290 K = J,N,MW
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
 280          DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) THEN
                IF (FAC(K) .LT. FACMAX) THEN
                  FAC(K) = MIN(100.D0*FAC(K), FACMAX)
                  GO TO 280
                ELSE
                  DY = YS
                END IF
              END IF
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              DFDY(MW,K) = Y(K)
 290          Y(K) = Y(K) + DY
            CALL F (N, T, Y, SAVE1)
            IF (N .EQ. 0) THEN
              JSTATE = 6
              RETURN
            END IF
            DO 330 K = J,N,MW
              Y(K) = DFDY(MW,K)
              YS = MAX(ABS(YWT(K)), ABS(Y(K)))
              DY = FAC(K)*YS
              IF (DY .EQ. 0.D0) DY = YS
              IF (NQ .EQ. 1) THEN
                DY = SIGN(DY, SAVE2(K))
              ELSE
                DY = SIGN(DY, YH(K,3))
              END IF
              DY = (Y(K) + DY) - Y(K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
C                                                                 Step 1
              IMAX = MAX(1, K - MU)
              DIFF = ABS(SAVE2(IMAX) - SAVE1(IMAX))
              I1 = IMAX
              I2 = MIN(K + ML, N)
              DO 310 I = I1+1,I2
                IF (ABS(SAVE2(I) - SAVE1(I)) .GT. DIFF) THEN
                  IMAX = I
                  DIFF = ABS(SAVE2(I) - SAVE1(I))
                END IF
 310            CONTINUE
C                                                                 Step 2
              IF (MIN(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX))) .GT.0.D0) THEN
                SCALE = MAX(ABS(SAVE2(IMAX)), ABS(SAVE1(IMAX)))
C                                                                 Step 3
                IF (DIFF .GT. BU*SCALE) THEN
                  FAC(K) = MAX(FACMIN, FAC(K)*.1D0)
                ELSE IF (BR*SCALE .LE.DIFF .AND. DIFF .LE.BL*SCALE) THEN
                  FAC(K) = MIN(FAC(K)*10.D0, FACMAX)
C                                                                 Step 4
                ELSE IF (DIFF .LT. BR*SCALE) THEN
                  FAC(K) = MIN(BP*FAC(K), FACMAX)
                END IF
              END IF
 330          CONTINUE
 340        CONTINUE
          NFE = NFE + J2
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.D0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 345 I = I1,I2
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
          BND = 0.D0
          IF (DFDYMX .NE. 0.D0) THEN
            DO 350 J = 1,N
              I1 = MAX(ML+1, MW+1-J)
              I2 = MIN(MW+N-J, MW+ML)
              DO 350 I = I1,I2
 350            BND = BND + (DFDY(I,J)/DFDYMX)**2
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          IF (N .EQ. 0) THEN
            JSTATE = 9
            RETURN
          END IF
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL DGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        IF (N .EQ. 0) THEN
          JSTATE = 10
          RETURN
        END IF
      END IF
      END
      SUBROUTINE DDRIV1 (N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDRIV1
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  880229   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV1 is to solve N (200 or fewer)
C            ordinary differential equations of the form
C            dY(I)/dT = F(Y(I),T), given the initial conditions
C            Y(I) = YI.  DDRIV1 uses double precision arithmetic.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C   Version 88.1
C
C  I.  CHOOSING THE CORRECT ROUTINE  ...................................
C
C     SDRIV
C     DDRIV
C     CDRIV
C           These are the generic names for three packages for solving
C           initial value problems for ordinary differential equations.
C           SDRIV uses single precision arithmetic.  DDRIV uses double
C           precision arithmetic.  CDRIV allows complex-valued
C           differential equations, integrated with respect to a single,
C           real precision, independent variable.
C
C    As an aid in selecting the proper program, the following is a
C    discussion of the important options or restrictions associated with
C    each program:
C
C      A. DDRIV1 should be tried first for those routine problems with
C         no more than 200 differential equations.  Internally this
C         routine has two important technical defaults:
C           1. Numerical approximation of the Jacobian matrix of the
C              right hand side is used.
C           2. The stiff solver option is used.
C         Most users of DDRIV1 should not have to concern themselves
C         with these details.
C
C      B. DDRIV2 should be considered for those problems for which
C         DDRIV1 is inadequate (SDRIV2 has no explicit restriction on
C         the number of differential equations.)  For example, DDRIV1
C         may have difficulty with problems having zero initial
C         conditions and zero derivatives.  In this case DDRIV2, with an
C         appropriate value of the parameter EWT, should perform more
C         efficiently.  DDRIV2 provides three important additional
C         options:
C           1. The nonstiff equation solver (as well as the stiff
C              solver) is available.
C           2. The root-finding option is available.
C           3. The program can dynamically select either the non-stiff
C              or the stiff methods.
C         Internally this routine also defaults to the numerical
C         approximation of the Jacobian matrix of the right hand side.
C
C      C. DDRIV3 is the most flexible, and hence the most complex, of
C         the programs.  Its important additional features include:
C           1. The ability to exploit band structure in the Jacobian
C              matrix.
C           2. The ability to solve some implicit differential
C              equations, i.e., those having the form:
C                   A(Y,T)*dY/dT = F(Y,T).
C           3. The option of integrating in the one step mode.
C           4. The option of allowing the user to provide a routine
C              which computes the analytic Jacobian matrix of the right
C              hand side.
C           5. The option of allowing the user to provide a routine
C              which does all the matrix algebra associated with
C              corrections to the solution components.
C
C  II.  ABSTRACT  ......................................................
C
C    The function of DDRIV1 is to solve N (200 or fewer) ordinary
C    differential equations of the form dY(I)/dT = F(Y(I),T), given the
C    initial conditions Y(I) = YI.  DDRIV1 is to be called once for each
C    output point.
C
C  III.  PARAMETERS  ...................................................
C
C       (REMEMBER--To run DDRIV1 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
c
C    The user should use parameter names in the call sequence of DDRIV1
C    for those quantities whose value may be altered by DDRIV1.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations, N .LE. 200
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routine F.  Thus
C             parameters required by F can be stored in this array in
C             components N+1 and above.  (Note: Changes by the user to
C             the first N components of this array will take effect only
C             after a restart, i.e., after setting MSTATE to +1(-1).)
C
C    TOUT   = (Input) The point at which the solution is desired.
C
C    MSTATE = An integer describing the status of integration.  The user
C             must initialize MSTATE to +1 or -1.  If MSTATE is
C             positive, the routine will integrate past TOUT and
C             interpolate the solution.  This is the most efficient
C             mode.  If MSTATE is negative, the routine will adjust its
C             internal step to reach TOUT exactly (useful if a
C             singularity exists beyond TOUT.)  The meaning of the
C             magnitude of MSTATE:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of MSTATE should be tested by the
C                  user.  Unless DDRIV1 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, DDRIV1 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV1
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV1
C                  again.
C               5  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.  See description of F in Section IV.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  On output, the adjusted relative accuracy if
C             the input value was too small.  The value of EPS should be
C             set as large as is reasonable, because the amount of work
C             done by DDRIV1 increases as EPS decreases.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The length of WORK should be at least N*N + 11*N + 225
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to DDRIV1.
C
C***LONG DESCRIPTION
C
C  IV.  USAGE  .........................................................
C
C                   PROGRAM SAMPLE
C                   DOUBLE PRECISION ALFA, EPS, T, TOUT
CC                                          N is the number of equations
C                   PARAMETER(ALFA = 1.D0, N = 3,
C                  8          LENW = N*N + 11*N + 225)
C                   DOUBLE PRECISION WORK(LENW), Y(N+1)
CC                                                         Initial point
C                   T = 0.00001D0
CC                                                Set initial conditions
C                   Y(1) = 10.D0
C                   Y(2) = 0.D0
C                   Y(3) = 10.D0
CC                                                        Pass parameter
C                   Y(4) = ALFA
C                   TOUT = T
C                   MSTATE = 1
C                   EPS = .001D0
C              10   CALL DDRIV1 (N, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C                   IF (MSTATE .GT. 2) STOP
C                   WRITE(*, '(4E12.3)') TOUT, (Y(I), I=1,3)
C                   TOUT = 10.D0*TOUT
C                   IF (TOUT .LT. 50.D0) GO TO 10
C                   END
C
C    The user must write a subroutine called F to evaluate the right
C    hand side of the differential equations.  It is of the form:
C
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION ALFA, T, Y(*), YDOT(*)
C                   ALFA = Y(N+1)
C                   YDOT(1) = 1.D0 + ALFA*(Y(2) - Y(1)) - Y(1)*Y(3)
C                   YDOT(2) = ALFA*(Y(1) - Y(2)) - Y(2)*Y(3)
C                   YDOT(3) = 1.D0 - Y(3)*(Y(1) + Y(2))
C                   END
C
C    This computes YDOT = F(Y,T), the right hand side of the
C    differential equations.  Here Y is a vector of length at least N.
C    The actual length of Y is determined by the user's declaration in
C    the program which calls DDRIV1.  Thus the dimensioning of Y in F,
C    while required by FORTRAN convention, does not actually allocate
C    any storage.  When this subroutine is called, the first N
C    components of Y are intermediate approximations to the solution
C    components.  The user should not alter these values.  Here YDOT is
C    a vector of length N.  The user should only compute YDOT(I) for I
C    from 1 to N.  Normally a return from F passes control back to
C    DDRIV1.  However, if the user would like to abort the calculation,
C    i.e., return control to the program which calls DDRIV1, he should
C    set N to zero.  DDRIV1 will signal this by returning a value of
C    MSTATE equal to +5(-5).  Altering the value of N in F has no effect
C    on the value of N in the call sequence of DDRIV1.
C
C  V.  OTHER COMMUNICATION TO THE USER  ................................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERROR.  That program will
C       terminate the user's run if it detects a probable problem setup
C       error, e.g., insufficient storage allocated by the user for the
C       WORK array.  For further information see section III-A of the
C       writeup for DDRIV3.
C
C    B. The number of evaluations of the right hand side can be found
C       in the WORK array in the location determined by:
C            LENW - (N + 21) + 4
C
C  VI.  REMARKS  .......................................................
C
C    For other information, see section IV of the writeup for DDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,XERROR
C***END PROLOGUE  DDRIV1
      EXTERNAL F
      DOUBLE PRECISION EPS, EWT, HMAX, T, TOUT, WORK(*), Y(*)
      PARAMETER(MXN = 200, IDLIW = 21)
      INTEGER IWORK(IDLIW+MXN)
      CHARACTER MSG*103
      PARAMETER(NROOT = 0, EWT = 1.D0, IERROR = 2, MINT = 2, MITER = 2,
     8          IMPL = 0, MXORD = 5, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  DDRIV1
      IF (N .GT. MXN) THEN
        WRITE(MSG, '(''DDRIV115FE Illegal input.  The number of '',
     8  ''equations,'', I8, '', is greater than the maximum allowed.'')
     8  ') N
        CALL XERROR(MSG(1:97), 97, 15, 2)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = 2.D0*ABS(TOUT - T)
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 10*N + 204)) THEN
        LNWCHK = N*N + 10*N + 204 + LENIW
        WRITE(MSG, '(''DDRIV116FE Insufficient storage allocated for '',
     8  ''the work array.  The required storage is at least'', I8)')
     8  LNWCHK
        CALL XERROR(MSG(1:103), 103, 16, 2)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
          II = I + LENWCM
 20       IWORK(I) = INT(WORK(II))
      END IF
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWT,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, F, F, NDE, MXSTEP, F, F)
      DO 40 I = 1,LENIW
        II = LENWCM + I
 40     WORK(II) = DBLE(IWORK(I))
      IF (NSTATE .LE. 4) THEN
        MSTATE = SIGN(NSTATE, MSTATE)
      ELSE IF (NSTATE .EQ. 6) THEN
        MSTATE = SIGN(5, MSTATE)
      END IF
      END
      SUBROUTINE DDRIV2 (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     8   LENW,IWORK,LENIW,G)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDRIV2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  871105   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV2 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  DDRIV2 uses double precision arithmetic.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C  I.  ABSTRACT  .......................................................
C
C    The function of DDRIV2 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.
C    DDRIV2 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C       (REMEMBER--To run DDRIV2 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C
C    The user should use parameter names in the call sequence of DDRIV2
C    for those quantities whose value may be altered by DDRIV2.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of differential equations.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F and
C             G.  Thus parameters required by F and G can be stored in
C             this array in components N+1 and above.  (Note: Changes
C             by the user to the first N components of this array will
C             take effect only after a restart, i.e., after setting
C             MSTATE to +1(-1).)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV2.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  DDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV2, he should set N to zero.
C             DDRIV2 will signal this by returning a value of MSTATE
C             equal to +6(-6).  Altering the value of N in F has no
C             effect on the value of N in the call sequence of DDRIV2.
C
C    TOUT   = (Input) The point at which the solution is desired.
C
C    MSTATE = An integer describing the status of integration.  The user
C             must initialize MSTATE to +1 or -1.  If MSTATE is
C             positive, the routine will integrate past TOUT and
C             interpolate the solution.  This is the most efficient
C             mode.  If MSTATE is negative, the routine will adjust its
C             internal step to reach TOUT exactly (useful if a
C             singularity exists beyond TOUT.)  The meaning of the
C             magnitude of MSTATE:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of MSTATE should be tested by the
C                  user.  Unless DDRIV2 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, DDRIV2 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV2
C                  again.  Other than an error in problem setup, the
C                  most likely cause for this condition is trying to
C                  integrate a stiff set of equations with the non-stiff
C                  integrator option. (See description of MINT below.)
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV2
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV2 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV2 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV2 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The index
C             of the equation whose root is being reported is stored in
C             the sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV2 increases as
C             EPS decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest physically
C             meaningful value for the solution.  This is used inter-
C             nally to compute an array YWT(I) = MAX(ABS(Y(I)), EWT).
C             One step error estimates divided by YWT(I) are kept less
C             than EPS.  Setting EWT to zero provides pure relative
C             error control.  However, setting EWT smaller than
C             necessary can adversely affect the running time.
C
C    MINT   = (Input) The integration method flag.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.
C             MINT may not be changed without restarting, i.e., setting
C             the magnitude of MSTATE to 1.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The length of WORK should be at least
C               16*N + 2*NROOT + 204         if MINT is 1, or
C               N*N + 10*N + 2*NROOT + 204   if MINT is 2, or
C               N*N + 17*N + 2*NROOT + 204   if MINT is 3,
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to DDRIV2.
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               21      if MINT is 1, or
C               N+21    if MINT is 2 or 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV2.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to
C             obtain the value of each of the NROOT equations for which
C             a root is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV2.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  DDRIV2.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV2, he should set N to zero.
C             DDRIV2 will signal this by returning a value of MSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of DDRIV2.
C
C***LONG DESCRIPTION
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERROR.  That program will
C       terminate the user's run if it detects a probable problem setup
C       error, e.g., insufficient storage allocated by the user for the
C       WORK array.  Messages are written on the standard error message
C       file.  At installations which have this error handling package
C       the user should determine the standard error handling file from
C       the local documentation.  Otherwise the short but serviceable
C       routine, XERROR, available with this package, can be used.  That
C       program writes on logical unit 6 to transmit messages.  A
C       complete description of XERROR is given in the Sandia
C       Laboratories report SAND78-1189 by R. E. Jones.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken since last initialization.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. On any return from DDRIV2 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    B. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV2.
C
C    C. When the routine G is not required, difficulties associated with
C       an unsatisfied external can be avoided by using the name of the
C       routine which calculates the right hand side of the differential
C       equations in place of G in the call sequence of DDRIV2.
C
C  V.  USAGE  ..........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               PARAMETER(MINT = 1, NROOT = 0, N = ...,
C              8          LENW = 16*N + 2*NROOT + 204, LENIW = 21)
C                                           N is the number of equations
C               DOUBLE PRECISION EPS, EWT, T, TOUT, WORK(LENW), Y(N)
C               INTEGER IWORK(LENIW)
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C               T = 0.                           Initial point
C               DO 10 I = 1,N
C          10     Y(I) = ...                     Set initial conditions
C               TOUT = T
C               EWT = ...
C               MSTATE = 1
C               EPS = ...
C          20   CALL DDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, F)
C                                          Last argument is not the same
C                                          as F if rootfinding is used.
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,XERROR
C***END PROLOGUE  DDRIV2
      EXTERNAL F, G
      DOUBLE PRECISION EPS, EWT, EWTCOM(1), G, HMAX, T, TOUT,
     8     WORK(*), Y(*)
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  DDRIV2
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(MSG, '(''DDRIV21FE Illegal input.  Improper value for '',
     8  ''the integration method flag,'', I8)') MINT
        CALL XERROR(MSG(1:81), 81, 21, 2)
        RETURN
      END IF
      IF (MSTATE .GE. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      EWTCOM(1) = EWT
      IF (EWT .NE. 0.D0) THEN
        IERROR = 3
      ELSE
        IERROR = 2
      END IF
      IF (MINT .EQ. 1) THEN
        MITER = 0
        MXORD = 12
      ELSE IF (MINT .EQ. 2) THEN
        MITER = 2
        MXORD = 5
      ELSE IF (MINT .EQ. 3) THEN
        MITER = 2
        MXORD = 12
      END IF
      HMAX = 2.D0*ABS(TOUT - T)
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, F, F, NDE, MXSTEP, G, F)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE DDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G,USERS)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDRIV3
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  880308   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             DOUBLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of DDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  DDRIV3
C            uses double precision arithmetic.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of DDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, DDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    DDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C       (REMEMBER--To run DDRIV3 correctly in double precision, ALL
C       non-integer arguments in the call sequence, including
C       arrays, MUST be declared double precision.)
C
C    The user should use parameter names in the call sequence of DDRIV3
C    for those quantities whose value may be altered by DDRIV3.  The
C    parameters in the call sequence are:
C
C    N      = (Input) The number of dependent functions whose solution
C             is desired.  N must not be altered during a problem.
C
C    T      = The independent variable.  On input for the first call, T
C             is the initial point.  On output, T is the point at which
C             the solution is given.
C
C    Y      = The vector of dependent variables.  Y is used as input on
C             the first call, to set the initial values.  On output, Y
C             is the computed solution vector.  This array Y is passed
C             in the call sequence of the user-provided routines F,
C             JACOBN, FA, USERS, and G.  Thus parameters required by
C             those routines can be stored in this array in components
C             N+1 and above.  (Note: Changes by the user to the first
C             N components of this array will take effect only after a
C             restart, i.e., after setting NSTATE to 1 .)
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C             This computes YDOT = F(Y,T), the right hand side of the
C             differential equations.  Here Y is a vector of length at
C             least N.  The actual length of Y is determined by the
C             user's declaration in the program which calls DDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.  Normally a return from F passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to 6 .  Altering the value of N in F has no effect
C             on the value of N in the call sequence of DDRIV3.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, DDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling DDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling DDRIV3 again.
C               6  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE F.
C               7  (Output)(Unsuccessful) N has been set to zero in
C                  FUNCTION G.  See description of G below.
C               8  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE JACOBN.  See description of JACOBN below.
C               9  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE FA.  See description of FA below.
C              10  (Output)(Unsuccessful) N has been set to zero in
C                  SUBROUTINE USERS.  See description of USERS below.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means DDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means DDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means DDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  DDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  DDRIV3 can find at most one root per
C             equation per internal integration step, and will then
C             return the solution either at TOUT or at a root, whichever
C             occurs first in the direction of integration.  The index
C             of the equation whose root is being reported is stored in
C             the sixth element of IWORK.
C             NOTE: NROOT is never altered by this program.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  EPS = 0 is allowed.  On output, the adjusted
C             relative accuracy if the input value was too small.  The
C             value of EPS should be set as large as is reasonable,
C             because the amount of work done by DDRIV3 increases as EPS
C             decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
C             physically meaningful value for the solution.  (Array,
C             possibly of length one.  See following description of
C             IERROR.)  Setting EWT smaller than necessary can adversely
C             affect the running time.
C
C    IERROR = (Input) Error control indicator.  A value of 3 is
C             suggested for most problems.  Other choices and detailed
C             explanations of EWT and IERROR are given below for those
C             who may need extra flexibility.
C
C             These last three input quantities EPS, EWT and IERROR
C             control the accuracy of the computed solution.  EWT and
C             IERROR are used internally to compute an array YWT.  One
C             step error estimates divided by YWT(I) are kept less than
C             EPS in root mean square norm.
C                 IERROR (Set by the user) =
C                 1  Means YWT(I) = 1. (Absolute error control)
C                                   EWT is ignored.
C                 2  Means YWT(I) = ABS(Y(I)),  (Relative error control)
C                                   EWT is ignored.
C                 3  Means YWT(I) = MAX(ABS(Y(I)), EWT(1)).
C                 4  Means YWT(I) = MAX(ABS(Y(I)), EWT(I)).
C                    This choice is useful when the solution components
C                    have differing scales.
C                 5  Means YWT(I) = EWT(I).
C             If IERROR is 3, EWT need only be dimensioned one.
C             If IERROR is 4 or 5, the user must dimension EWT at least
C             N, and set its values.
C
C    MINT   = (Input) The integration method indicator.
C               MINT = 1  Means the Adams methods, and is used for
C                         non-stiff problems.
C               MINT = 2  Means the stiff methods of Gear (i.e., the
C                         backward differentiation formulas), and is
C                         used for stiff problems.
C               MINT = 3  Means the program dynamically selects the
C                         Adams methods when the problem is non-stiff
C                         and the Gear methods when the problem is
C                         stiff.  When using the Adams methods, the
C                         program uses a value of MITER=0; when using
C                         the Gear methods, the program uses the value
C                         of MITER provided by the user.  Only a value
C                         of IMPL = 0 and a value of MITER = 1, 2, 4, or
C                         5 is allowed for this option.  The user may
C                         not alter the value of MINT or MITER without
C                         restarting, i.e., setting NSTATE to 1.
C
C    MITER  = (Input) The iteration method indicator.
C               MITER = 0  Means functional iteration.  This value is
C                          suggested for non-stiff problems.
C               MITER = 1  Means chord method with analytic Jacobian.
C                          In this case, the user supplies subroutine
C                          JACOBN (see description below).
C               MITER = 2  Means chord method with Jacobian calculated
C                          internally by finite differences.
C               MITER = 3  Means chord method with corrections computed
C                          by the user-written routine USERS (see
C                          description of USERS below.)  This option
C                          allows all matrix algebra and storage
C                          decisions to be made by the user.  When using
C                          a value of MITER = 3, the subroutine FA is
C                          not required, even if IMPL is not 0.  For
C                          further information on using this option, see
C                          section IV-E below.
C               MITER = 4  Means the same as MITER = 1 but the A and
C                          Jacobian matrices are assumed to be banded.
C               MITER = 5  Means the same as MITER = 2 but the A and
C                          Jacobian matrices are assumed to be banded.
C
C    IMPL   = (Input) The implicit method indicator.
C               IMPL = 0 Means solving dY(I)/dT = F(Y(I),T).
C               IMPL = 1 Means solving A*dY(I)/dT = F(Y(I),T),
C                        non-singular A (see description of FA below.)
C                        Only MINT = 1 or 2, and MITER = 1, 2, 3, 4, or
C                        5 are allowed for this option.
C               IMPL = 2 Means solving certain systems of hybrid
C                        differential/algebraic equations (see
C                        description of FA below.)  Only MINT = 2 and
C                        MITER = 1, 2, 3, 4, or 5, are allowed for this
C                        option.
C               The value of IMPL must not be changed during a problem.
C
C    ML     = (Input) The lower half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(R-C) for nonzero
C             A(R,C).)
C
C    MU     = (Input) The upper half-bandwidth in the case of a banded
C             A or Jacobian matrix.  (I.e., maximum(C-R).)
C
C    MXORD  = (Input) The maximum order desired. This is .LE. 12 for
C             the Adams methods and .LE. 5 for the Gear methods.  Normal
C             value is 12 and 5, respectively.  If MINT is 3, the
C             maximum order used will be MIN(MXORD, 12) when using the
C             Adams methods, and MIN(MXORD, 5) when using the Gear
C             methods.  MXORD must not be altered during a problem.
C
C    HMAX   = (Input) The maximum magnitude of the step size that will
C             be used for the problem.  This is useful for ensuring that
C             important details are not missed.  If this is not the
C             case, a large value, such as the interval length, is
C             suggested.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW double precision words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       DOUBLE PRECISION WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             DDRIV3.
C
C      IMPL =   0                   1                   2
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N +       Not allowed         Not allowed
C              2*NROOT + 204
C
C         1,2  N*N+(MXORD+5)*N     2*N*N+(MXORD+5)*N   N*N+(MXORD+6)*N
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204
C
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +
C              (MXORD+6)*N +       (MXORD+7)*N +       (MXORD+7)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C              ---------------------------------------------------------
C
C    IWORK
C    LENIW  = (Input)
C             IWORK is an integer array of length LENIW used internally
C             for temporary storage.  The user must allocate space for
C             this array in the calling program by a statement such as
C                       INTEGER IWORK(...)
C             The length of IWORK should be at least
C               21      if MITER is 0 or 3, or
C               N+21    if MITER is 1, 2, 4, or 5, or MINT is 3,
C             and LENIW should be set to the value used.  The contents
C             of IWORK should not be disturbed between calls to DDRIV3.
C
C    JACOBN = A subroutine supplied by the user, if MITER is 1 or 4.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  Given a system of N
C             differential equations, it is meaningful to speak about
C             the partial derivative of the I-th right hand side with
C             respect to the J-th dependent variable.  In general there
C             are N*N such quantities.  Often however the equations can
C             be ordered so that the I-th differential equation only
C             involves dependent variables with index near I, e.g., I+1,
C             I-2.  Such a system is called banded.  If, for all I, the
C             I-th equation depends on at most the variables
C               Y(I-ML), Y(I-ML+1), ... , Y(I), Y(I+1), ... , Y(I+MU)
C             then we call ML+MU+1 the bandwith of the system.  In a
C             banded system many of the partial derivatives above are
C             automatically zero.  For the cases MITER = 1, 2, 4, and 5,
C             some of these partials are needed.  For the cases
C             MITER = 2 and 5 the necessary derivatives are
C             approximated numerically by DDRIV3, and we only ask the
C             user to tell DDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   DOUBLE PRECISION Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             JACOBN, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  If the system is not
C             banded (MITER=1), the partials of the I-th equation with
C             respect to the J-th dependent function are to be stored in
C             DFDY(I,J).  Thus partials of the I-th equation are stored
C             in the I-th row of DFDY.  If the system is banded
C             (MITER=4), then the partials of the I-th equation with
C             respect to Y(J) are to be stored in DFDY(K,J), where
C             K=I-J+MU+1 .  Normally a return from JACOBN passes control
C             back to DDRIV3.  However, if the user would like to abort
C             the calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +8(-8).
C             Altering the value of N in JACOBN has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C    FA     = A subroutine supplied by the user if IMPL is 1 or 2, and
C             MITER is not 3.  If so, the name must be declared EXTERNAL
C             in the user's calling program.  This subroutine computes
C             the array A, where A*dY(I)/dT = F(Y(I),T).
C             There are two cases:
C
C               IMPL=1.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(MATDIM,*)
C                     .
C                     .
C                     Calculate ALL values of A
C                     .
C                     .
C                   END (Sample)
C               In this case A is assumed to be a nonsingular matrix,
C               with the same structure as DFDY (see JACOBN description
C               above).  Programming considerations prevent complete
C               generality.  If MITER is 1 or 2, A is assumed to be full
C               and the user must compute and store all values of
C               A(I,J), I,J=1, ... ,N.  If MITER is 4 or 5, A is assumed
C               to be banded with lower and upper half bandwidth ML and
C               MU.  The left hand side of the I-th equation is a linear
C               combination of dY(I-ML)/dT, dY(I-ML+1)/dT, ... ,
C               dY(I)/dT, ... , dY(I+MU-1)/dT, dY(I+MU)/dT.  Thus in the
C               I-th equation, the coefficient of dY(J)/dT is to be
C               stored in A(K,J), where K=I-J+MU+1.
C               NOTE: The array A will be altered between calls to FA.
C
C               IMPL=2.
C               Subroutine FA is of the form:
C                   SUBROUTINE FA (N, T, Y, A, MATDIM, ML, MU, NDE)
C                   DOUBLE PRECISION Y(*), A(*)
C                     .
C                     .
C                     Calculate non-zero values of A(1),...,A(NDE)
C                     .
C                     .
C                   END (Sample)
C               In this case it is assumed that the system is ordered by
C               the user so that the differential equations appear
C               first, and the algebraic equations appear last.  The
C               algebraic equations must be written in the form:
C               0 = F(Y(I),T).  When using this option it is up to the
C               user to provide initial values for the Y(I) that satisfy
C               the algebraic equations as well as possible.  It is
C               further assumed that A is a vector of length NDE.  All
C               of the components of A, which may depend on T, Y(I),
C               etc., must be set by the user to non-zero values.
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.  Normally a return from FA passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +9(-9).
C             Altering the value of N in FA has no effect on the value
C             of N in the call sequence of DDRIV3.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to DDRIV3.
C
C    G      = A double precision FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   DOUBLE PRECISION FUNCTION G (N, T, Y, IROOT)
C                   DOUBLE PRECISION Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls DDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.  Normally a return from G passes
C             control back to  DDRIV3.  However, if the user would like
C             to abort the calculation, i.e., return control to the
C             program which calls DDRIV3, he should set N to zero.
C             DDRIV3 will signal this by returning a value of NSTATE
C             equal to +7(-7).  In this case, the index of the equation
C             being evaluated is stored in the sixth element of IWORK.
C             Altering the value of N in G has no effect on the value of
C             N in the call sequence of DDRIV3.
C
C    USERS  = A subroutine supplied by the user, if MITER is 3.
C             If this is the case, the name must be declared EXTERNAL in
C             the user's calling program.  The routine USERS is called
C             by DDRIV3 when certain linear systems must be solved.  The
C             user may choose any method to form, store and solve these
C             systems in order to obtain the solution result that is
C             returned to DDRIV3.  In particular, this allows sparse
C             matrix methods to be used.  The call sequence for this
C             routine is:
C
C                SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL,
C               8                  IMPL, N, NDE, IFLAG)
C                DOUBLE PRECISION Y(*), YH(*), YWT(*), SAVE1(*),
C               8     SAVE2(*), T, H, EL
C
C             The input variable IFLAG indicates what action is to be
C             taken.  Subroutine USERS should perform the following
C             operations, depending on the value of IFLAG and IMPL.
C
C               IFLAG = 0
C                 IMPL = 0.  USERS is not called.
C                 IMPL = 1 or 2.  Solve the system A*X = SAVE2,
C                   returning the result in SAVE2.  The array SAVE1 can
C                   be used as a work array.
C
C               IFLAG = 1
C                 IMPL = 0.  Compute, decompose and store the matrix
C                   (I - H*EL*J), where I is the identity matrix and J
C                   is the Jacobian matrix of the right hand side.  The
C                   array SAVE1 can be used as a work array.
C                 IMPL = 1 or 2. Compute, decompose and store the matrix
C                   (A - H*EL*J).  The array SAVE1 can be used as a work
C                   array.
C
C               IFLAG = 2
C                 IMPL = 0.   Solve the system
C                     (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                   returning the result in SAVE2.
C                 IMPL = 1 or 2.  Solve the system
C                   (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                   returning the result in SAVE2.
C                 The array SAVE1 should not be altered.
C             Normally a return from USERS passes control back to
C             DDRIV3.  However, if the user would like to abort the
C             calculation, i.e., return control to the program which
C             calls DDRIV3, he should set N to zero.  DDRIV3 will signal
C             this by returning a value of NSTATE equal to +10(-10).
C             Altering the value of N in USERS has no effect on the
C             value of N in the call sequence of DDRIV3.
C
C***LONG DESCRIPTION
C
C  III.  OTHER COMMUNICATION TO THE USER  ..............................
C
C    A. The solver communicates to the user through the parameters
C       above.  In addition it writes diagnostic messages through the
C       standard error handling program XERROR.  That program will
C       terminate the user's run if it detects a probable problem setup
C       error, e.g., insufficient storage allocated by the user for the
C       WORK array.  Messages are written on the standard error message
C       file.  At installations which have this error handling package
C       the user should determine the standard error handling file from
C       the local documentation.  Otherwise the short but serviceable
C       routine, XERROR, available with this package, can be used.  That
C       program writes on logical unit 6 to transmit messages.  A
C       complete description of XERROR is given in the Sandia
C       Laboratories report SAND78-1189 by R. E. Jones.  Following is a
C       list of possible errors.  Unless otherwise noted, all messages
C       come from DDRIV3:
C
C        No.  Type         Message
C        ---  ----         -------
C         1   Fatal        From DDRIV2: The integration method flag has
C                          an illegal value.
C         2   Warning      The output point is inconsistent with the
C                          value of NTASK and T.
C         3   Warning      Number of steps to reach TOUT exceeds MXSTEP.
C         4   Recoverable  Requested accuracy is too stringent.
C         5   Warning      Step size is below the roundoff level.
C         6   Fatal        EPS is less than zero.
C         7   Fatal        N is not positive.
C         8   Fatal        Insufficient work space provided.
C         9   Fatal        Improper value for NSTATE, MINT, MITER and/or
C                          IMPL.
C        10   Fatal        The IWORK array is too small.
C        11   Fatal        The step size has gone to zero.
C        12   Fatal        Excessive amount of work.
C        13   Fatal        For IMPL=1 or 2, the matrix A is singular.
C        14   Fatal        MXORD is not positive.
C        15   Fatal        From DDRIV1: N is greater than 200.
C        16   Fatal        From DDRIV1: The WORK array is too small.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken since last initialization.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         DDNTP, DDZRO, DDSTP, DDNTL, DDPST, DDCOR, DDCST,
C         DDPSC, and DDSCL;
C         DGEFA, DGESL, DGBFA, DGBSL, and DNRM2 (from LINPACK)
C         D1MACH (from the Bell Laboratories Machine Constants Package)
C         XERROR (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from DDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV3.
C
C    D. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to DDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to DDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    E. As the price for complete control of matrix algebra, the DDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C              DOUBLE PRECISION DFDY(N,N), EPSJ, H, R, D1MACH,
C             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
C              UROUND = D1MACH(4)
C              EPSJ = SQRT(UROUND)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.D0) R = YWT(J)
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
C                IF (N .EQ. 0) RETURN
C                Y(J) = YJ
C                DO 20 I = I1,I2
C         20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
C         30     CONTINUE
C       Many problems give rise to structured sparse Jacobians, e.g.,
C       block banded.  It is possible to approximate them with fewer
C       function evaluations than the above procedure uses; see Curtis,
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
C       pp. 117-119.
C
C    F. When any of the routines JACOBN, FA, G, or USERS, is not
C       required, difficulties associated with unsatisfied externals can
C       be avoided by using the name of the routine which calculates the
C       right hand side of the differential equations in place of the
C       corresponding name in the call sequence of DDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDSTP,DDNTP,DDZRO,DGEFA,DGESL,DGBFA,DGBSL,DNRM2,
C                    D1MACH,XERROR
C***END PROLOGUE  DDRIV3
      EXTERNAL F, JACOBN, FA, G, USERS
      DOUBLE PRECISION AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,
     8     NROUND, RE, D1MACH, SIZE, DNRM2, SUM, T, TLAST, TOUT, TROOT,
     8     UROUND, WORK(*), Y(*)
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.D0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,
     8          INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20, INDPRT = 21,
     8          INDPVT = 22)
C***FIRST EXECUTABLE STATEMENT  DDRIV3
      NPAR = N
      UROUND = D1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = D1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.D0) THEN
        WRITE(MSG, '(''DDRIV36FE Illegal input.  EPS,'', D16.8,
     8  '', is negative.'')') EPS
        CALL XERROR(MSG(1:60), 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''DDRIV37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
        CALL XERROR(MSG(1:72), 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''DDRIV314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
        CALL XERROR(MSG(1:67), 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1) .OR.
     8  (NSTATE .LT. 1 .OR. NSTATE .GT. 10)) THEN
        WRITE(MSG, '(''DDRIV39FE Illegal input.  Improper value for '',
     8  ''NSTATE(MSTATE), MINT, MITER or IMPL.'')')
        CALL XERROR(MSG(1:81), 81, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(MSG, '(''DDRIV310FE Illegal input.  Insufficient '',
     8  ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LIWCHK
        CALL XERROR(MSG(1:164), 164, 10, 2)
        RETURN
      END IF
C                                                Allocate the WORK array
C                                         IYH is the index of YH in WORK
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
C                                             IDFDY is the index of DFDY
C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3)  THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2)  THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5)  THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
C                                           ITROOT is the index of TROOT
      IFAC = ITROOT + NROOT
C                                               IFAC is the index of FAC
      IF (MITER .EQ. 2 .OR. MITER .EQ. 5 .OR. MINT .EQ. 3) THEN
        IA = IFAC + N
      ELSE
        IA = IFAC
      END IF
C                                                   IA is the index of A
      IF (IMPL .EQ. 0 .OR. MITER .EQ. 3) THEN
        LENCHK = IA - 1
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
        LENCHK = IA - 1 + N*N
      ELSE IF (IMPL .EQ. 1 .AND. (MITER .EQ. 4 .OR. MITER .EQ. 5)) THEN
        LENCHK = IA - 1 + (2*ML + MU + 1)*N
      ELSE IF (IMPL .EQ. 2 .AND. MITER .NE. 3) THEN
        LENCHK = IA - 1 + N
      END IF
      IF (LENW .LT. LENCHK) THEN
        WRITE(MSG, '(''DDRIV38FE Illegal input.  Insufficient '',
     8  ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LENCHK
        CALL XERROR(MSG(1:162), 162, 8, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        MATDIM = 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        MATDIM = N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MATDIM = 2*ML + MU + 1
      END IF
      IF (IMPL .EQ. 0 .OR. IMPL .EQ. 1) THEN
        NDECOM = N
      ELSE IF (IMPL .EQ. 2) THEN
        NDECOM = NDE
      END IF
      IF (NSTATE .EQ. 1) THEN
C                                                  Initialize parameters
        IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
          IWORK(IMXORD) = MIN(MXORD, 12)
        ELSE IF (MINT .EQ. 2) THEN
          IWORK(IMXORD) = MIN(MXORD, 5)
        END IF
        IWORK(IMXRDS) = MXORD
        IF (MINT .EQ. 1 .OR. MINT .EQ. 2) THEN
          IWORK(IMNT) = MINT
          IWORK(IMTR) = MITER
          IWORK(IMNTLD) = MINT
          IWORK(IMTRLD) = MITER
        ELSE IF (MINT .EQ. 3) THEN
          IWORK(IMNT) = 1
          IWORK(IMTR) = 0
          IWORK(IMNTLD) = IWORK(IMNT)
          IWORK(IMTRLD) = IWORK(IMTR)
          IWORK(IMTRSV) = MITER
        END IF
        WORK(IHMAX) = HMAX
        H = (TOUT - T)*(1.D0 - 4.D0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.D0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.D0
        WORK(IHUSED) =0.D0
        WORK(IAVGRD) = 0.D0
        IWORK(INDMXR) = 0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        IWORK(INROOT) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
        IF (T .EQ. TOUT) RETURN
        GO TO 180
      END IF
C                                             On a continuation, check
C                                             that output points have
C                                             been or will be overtaken.
      IF (IWORK(ICNVRG) .EQ. 1) THEN
        CONVRG = .TRUE.
      ELSE
        CONVRG = .FALSE.
      END IF
      T = WORK(IT)
      H = WORK(IH)
      HSIGN = WORK(IHSIGN)
      IF (IWORK(IJTASK) .EQ. 0) GO TO 180
C
C                                   IWORK(IJROOT) flags unreported
C                                   roots, and is set to the value of
C                                   NTASK when a root was last selected.
C                                   It is set to zero when all roots
C                                   have been reported.  IWORK(INROOT)
C                                   contains the index and WORK(ITOUT)
C                                   contains the value of the root last
C                                   selected to be reported.
C                                   IWORK(INRTLD) contains the value of
C                                   NROOT and IWORK(INDTRT) contains
C                                   the value of ITROOT when the array
C                                   of roots was last calculated.
      IF (NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
              T = TROOT
              NSTATE = 5
              GO TO 580
            END IF
C                                         A root has just been reported.
C                                         Select the next root.
          ELSE
            TROOT = T
            IROOT = 0
            DO 50 I = 1,IWORK(INRTLD)
              JTROOT = IWORK(INDTRT) + I - 1
              IF (WORK(JTROOT)*HSIGN .LE. TROOT*HSIGN) THEN
C
C                                              Check for multiple roots.
C
                IF (WORK(JTROOT) .EQ. WORK(ITOUT) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                  GO TO 60
                END IF
                IF (WORK(JTROOT)*HSIGN .GT. WORK(ITOUT)*HSIGN) THEN
                  IROOT = I
                  TROOT = WORK(JTROOT)
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = TROOT
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
                  NSTATE = 5
                  T = TROOT
                  GO TO 580
                END IF
              END IF
            ELSE IF (NTASK .EQ. 2 .OR. NTASK .EQ. 3) THEN
C
C                                     If there are no more roots, or the
C                                     user has altered TOUT to be less
C                                     than a root, set IJROOT to zero.
C
              IF (IROOT .EQ. 0 .OR. (TOUT*HSIGN .LT. TROOT*HSIGN)) THEN
                IWORK(IJROOT) = 0
              ELSE
                CALL DDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
                NSTATE = 5
                T = TROOT
                GO TO 580
              END IF
            END IF
          END IF
        END IF
      END IF
C
      IF (NTASK .EQ. 1) THEN
        NSTATE = 2
        IF (T*HSIGN .GE. TOUT*HSIGN) THEN
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG(1:124), 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          NSTATE = 2
          GO TO 580
        END IF
C                                   Determine if TOUT has been overtaken
C
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          NSTATE = 2
          GO TO 560
        END IF
C                                             If there are no more roots
C                                             to report, report T.
        IF (NSTATE .EQ. 5) THEN
          NSTATE = 2
          GO TO 560
        END IF
        NSTATE = 2
C                                                       See if TOUT will
C                                                       be overtaken.
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''DDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', D16.8, '', was beyond TOUT,'', D16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG(1:124), 124, 2, 0)
          CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.D0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C                         Implement changes in MINT, MITER, and/or HMAX.
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. WORK(IHMAX)) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. WORK(IH)) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = H
        END IF
        WORK(IHMAX) = HMAX
      END IF
C
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
 190    Y(I) = WORK(JYH)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
 200     CONTINUE
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = 1.D0
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
 250      WORK(JYWT) = EWT(I)
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. 0.D0) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (NPAR, T, Y, WORK(ISAVE2))
          IF (NPAR .EQ. 0) THEN
            NSTATE = 6
            RETURN
          END IF
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),
     8                 WORK(ISAVE2), T, H, WORK(IEL), IMPL, NPAR,
     8                 NDECOM, IFLAG)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 10
              RETURN
            END IF
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (NPAR, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              IF (NPAR .EQ. 0) THEN
                NSTATE = 9
                RETURN
              END IF
              CALL DGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (NPAR, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            IF (NPAR .EQ. 0) THEN
              NSTATE = 9
              RETURN
            END IF
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF (WORK(JA) .EQ. 0.D0) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.D0) THEN
            WORK(JYWT) = ABS(Y(J))
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = ABS(H*WORK(JSAVE2))
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = ABS(WORK(JHYP))
            END IF
          END IF
          IF (WORK(JYWT) .EQ. 0.D0) WORK(JYWT) = UROUND
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
 380      WORK(JYWT) = MAX(EWT(1), ABS(Y(I)))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
 400      WORK(JYWT) = MAX(EWT(I), ABS(Y(I)))
      END IF
C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)
      SUM = DNRM2(N, WORK(ISAVE2), 1)/SQRT(DBLE(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.D0 + 10.D0*UROUND)
        WRITE(MSG, '(''DDRIV34REC At T,'', D16.8, '', the requested '',
     8  ''accuracy, EPS, was not obtainable with the machine '',
     8  ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(D16.8)') EPS
        CALL XERROR(MSG(1:152), 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(MSG, '(''DDRIV35WRN At T,'', D16.8, '', the step size,'',
     8  D16.8, '', is smaller than the roundoff level of T.  '')') T, H
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
     8  ''change in the right hand side of the differential '',
     8  ''equations.'')')
        CALL XERROR(MSG(1:205), 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''DDRIV33WRN At T,'', D16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', D16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG(1:103), 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND, USERS,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, FAC, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, NPAR,
     8           NDECOM, WORK(IYWT), UROUND, USERS,  WORK(IAVGH),
     8           WORK(IAVGRD), WORK(IH), WORK(IHUSED), IWORK(IJTASK),
     8           IWORK(IMNTLD), IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IFAC),
     8            WORK(IHOLD), IWORK(INDPVT), JSTATE, IWORK(INQ),
     8            IWORK(INWAIT), WORK(IRC), WORK(IRMAX), WORK(ISAVE1),
     8            WORK(ISAVE2), WORK(ITQ), WORK(ITREND), MINT,
     8            IWORK(IMTRSV), IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (NPAR, T, Y, I)
          IF (NPAR .EQ. 0) THEN
            IWORK(INROOT) = I
            NSTATE = 7
            RETURN
          END IF
          IF (GLAST*WORK(JGNOW) .GT. 0.D0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.D0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.D0) THEN
                WORK(JTROOT) = T + H
              ELSE
                IF (ABS(WORK(IHUSED)) .GE. UROUND*ABS(T)) THEN
                  TLAST = T - WORK(IHUSED)
                  IROOT = I
                  TROOT = T
                  CALL DDZRO (AE, G, H, NPAR, IWORK(INQ), IROOT, RE, T,
     8                        WORK(IYH), UROUND,  TROOT, TLAST,
     8                        WORK(JGNOW), GLAST,  Y)
                  DO 480 J = 1,N
  480               Y(J) = WORK(IYH + J -1)
                  IF (NPAR .EQ. 0) THEN
                    IWORK(INROOT) = I
                    NSTATE = 7
                    RETURN
                  END IF
                  WORK(JTROOT) = TROOT
                ELSE
                  WORK(JTROOT) = T
                  IROOT = I
                END IF
              END IF
            END IF
          END IF
 500      CONTINUE
        IF (IROOT .EQ. 0) THEN
          IWORK(IJROOT) = 0
C                                                  Select the first root
        ELSE
          IWORK(IJROOT) = NTASK
          IWORK(INRTLD) = NROOT
          IWORK(INDTRT) = ITROOT
          TROOT = T + H
          DO 510 I = 1,NROOT
            JTROOT = ITROOT + I - 1
            IF (WORK(JTROOT)*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = WORK(JTROOT)
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = TROOT
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL DDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
            NSTATE = 5
            T = TROOT
            GO TO 580
          END IF
        END IF
      END IF
C                               Test for NTASK condition to be satisfied
      NSTATE = 2
      IF (NTASK .EQ. 1) THEN
        IF (T*HSIGN .LT. TOUT*HSIGN) GO TO 260
        CALL DDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
C                               TOUT is assumed to have been attained
C                               exactly if T is within twenty roundoff
C                               units of TOUT, relative to max(TOUT, T).
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.D0 - 4.D0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.D0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
          GO TO 260
        END IF
      END IF
C                                      All returns are made through this
C                                      section.  IMXERR is determined.
 560  DO 570 I = 1,N
        JYH = I + IYH - 1
 570    Y(I) = WORK(JYH)
 580  IF (CONVRG) THEN
        IWORK(ICNVRG) = 1
      ELSE
        IWORK(ICNVRG) = 0
      END IF
      IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.D0
      IMXERR = 1
      IWORK(INDMXR) = IMXERR
      DO  590 I = 1,N
C                                            SIZE = ABS(ERROR(I)/YWT(I))
        JYWT = I + IYWT - 1
        JERROR = I + ISAVE1 - 1
        SIZE = ABS(WORK(JERROR)/WORK(JYWT))
        IF (BIG .LT. SIZE) THEN
          BIG = SIZE
          IMXERR = I
          IWORK(INDMXR) = IMXERR
        END IF
 590    CONTINUE
      RETURN
C
 660  NSTATE = JSTATE
      RETURN
C                                        Fatal errors are processed here
C
 670  WRITE(MSG, '(''DDRIV311FE At T,'', D16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL XERROR(MSG(1:129), 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''DDRIV312FE At T,'', D16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG(1:165), 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''DDRIV313FE At T,'', D16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG(1:74), 74, 13, 2)
      RETURN
      END
      SUBROUTINE DDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDSCL
C***REFER TO  DDRIV3
C   This subroutine rescales the YH array whenever the step size
C   is changed.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850319   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDSCL
      DOUBLE PRECISION H, HMAX, RC, RH, RMAX, R1, YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDSCL
      IF (H .LT. 1.D0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.D0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      END
      SUBROUTINE DDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,
     8   MITER,ML,MU,N,NDE,YWT,UROUND,USERS,AVGH,AVGORD,H,HUSED,JTASK,
     8   MNTOLD,MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,FAC,
     8   HOLD,IPVT,JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,
     8   MTRSV,MXRDSV)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDSTP
C***REFER TO  DDRIV3
C  DDSTP performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Communication with DDSTP is done with the following variables:
C
C    YH      An N by MAXORD+1 array containing the dependent variables
C              and their scaled derivatives.  MAXORD, the maximum order
C              used, is currently 12 for the Adams methods and 5 for the
C              Gear methods.  YH(I,J+1) contains the J-th derivative of
C              Y(I), scaled by H**J/factorial(J).  Only Y(I),
C              1 .LE. I .LE. N, need be set by the calling program on
C              the first entry.  The YH array should not be altered by
C              the calling program.  When referencing YH as a
C              2-dimensional array, use a column length of N, as this is
C              the value used in DDSTP.
C    DFDY    A block of locations used for partial derivatives if MITER
C              is not 0.  If MITER is 1 or 2 its length must be at least
C              N*N.  If MITER is 4 or 5 its length must be at least
C              (2*ML+MU+1)*N.
C    YWT     An array of N locations used in convergence and error tests
C    SAVE1
C    SAVE2   Arrays of length N used for temporary storage.
C    IPVT    An integer array of length N used by the linear system
C              solvers for the storage of row interchange information.
C    A       A block of locations used to store the matrix A, when using
C              the implicit method.  If IMPL is 1, A is a MATDIM by N
C              array.  If MITER is 1 or 2 MATDIM is N, and if MITER is 4
C              or 5 MATDIM is 2*ML+MU+1.  If IMPL is 2 its length is N.
C    JTASK   An integer used on input.
C              It has the following values and meanings:
C                 .EQ. 0  Perform the first step.  This value enables
C                         the subroutine to initialize itself.
C                .GT. 0  Take a new step continuing from the last.
C                         Assumes the last step was successful and
C                         user has not changed any parameters.
C                 .LT. 0  Take a new step with a new value of H and/or
C                         MINT and/or MITER.
C    JSTATE  A completion code with the following meanings:
C                1  The step was successful.
C                2  A solution could not be obtained with H .NE. 0.
C                3  A solution was not obtained in MXTRY attempts.
C                4  For IMPL .NE. 0, the matrix A is singular.
C              On a return with JSTATE .GT. 1, the values of T and
C              the YH array are as of the beginning of the last
C              step, and H is the last step size attempted.
C***ROUTINES CALLED  DDNTL,DDPST,DDCOR,DDPSC,DDSCL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870810   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDSTP
      EXTERNAL F, JACOBN, FA, USERS
      DOUBLE PRECISION A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, FAC(*), H, HMAX, HN, HOLD, HS, HUSED,
     8     NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM,
     8     SAVE1(*), SAVE2(*), DNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD,
     8     UROUND, Y(*), YH(N,*), YWT(*), Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3D0, BIAS2 = 1.2D0, BIAS3 = 1.4D0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3D0, RMFAIL = 2.D0,
     8          RMNORM = 10.D0, TRSHLD = 1.D0)
      DATA IER /.FALSE./
C***FIRST EXECUTABLE STATEMENT  DDSTP
      NSV = N
      BND = 0.D0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
      END IF
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL DDPSC (1, N, NQ,  YH)
      EVALJC = ((ABS(RC - 1.D0) .GT. RCTEST) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      IF (N .EQ. 0) THEN
        JSTATE = 6
        GO TO 430
      END IF
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, USERS, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, FAC, IER, IPVT, SAVE1, ISWFLG,
     8              BND, JSTATE)
        IF (N .EQ. 0) GO TO 430
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.D0
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = 0.D0
C                      Up to MXITER corrector iterations are taken.
C                      Convergence is tested by requiring the r.m.s.
C                      norm of changes to be less than EPS.  The sum of
C                      the corrections is accumulated in the vector
C                      SAVE1(I).  It is approximately equal to the L-th
C                      derivative of Y multiplied by
C                      H**L/(factorial(L-1)*EL(L,NQ)), and is thus
C                      proportional to the actual errors to the lowest
C                      power of H present (H**L).  The YH array is not
C                      altered in the correction loop.  The norm of the
C                      iterate difference is stored in D.  If
C                      ITER .GT. 0, an estimate of the convergence rate
C                      constant is stored in TREND, and this is used in
C                      the convergence test.
C
 130  CALL DDCOR (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,
     8            MU, N, NDE, NQ, T, USERS, Y, YH, YWT,  EVALFA, SAVE1,
     8            SAVE2,  A, D, JSTATE)
        IF (N .EQ. 0) GO TO 430
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = DNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = DNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = DNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.D0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.D0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9D0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.D0*TREND, 1.D0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
        IF (N .EQ. 0) THEN
          JSTATE = 6
          GO TO 430
        END IF
        NFE = NFE + 1
        GO TO 130
      END IF
C                     The corrector iteration failed to converge in
C                     MXITER tries.  If partials are involved but are
C                     not up to date, they are reevaluated for the next
C                     try.  Otherwise the YH array is retracted to its
C                     values before prediction, and H is reduced, if
C                     possible.  If not, a no-convergence exit is taken.
      IF (CONVRG) THEN
        EVALJC = .TRUE.
        EVALFA = .FALSE.
        GO TO 110
      END IF
 160  T = TOLD
      CALL DDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3D0
      ELSE
        RH = .9D0*(EPS/CTEST)**(.2D0)
      END IF
      IF (RH*H .EQ. 0.D0) GO TO 400
      CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C                          The corrector has converged.  CONVRG is set
C                          to .TRUE. if partial derivatives were used,
C                          to indicate that they may need updating on
C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      DO 180 I = 1,NDE
 180    SAVE2(I) = SAVE1(I)/YWT(I)
      ETEST = DNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(DBLE(NDE)))
C
C                           The error test failed.  NFAIL keeps track of
C                           multiple failures.  Restore T and the YH
C                           array to their previous values, and prepare
C                           to try the step again.  Compute the optimum
C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL DDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.D0/(BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
          IF (NQ .GT. 1) THEN
            DO 190 I = 1,NDE
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
            RH1 = 1.D0/MAX(1.D0, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
            IF (RH2 .LT. RH1) THEN
              NQ = NQ - 1
              RC = RC*EL(1,NQ)/EL(1,NQ+1)
              RH = RH1
            ELSE
              RH = RH2
            END IF
          ELSE
            RH = RH2
          END IF
          NWAIT = NQ + 2
          IF (RH*H .EQ. 0.D0) GO TO 400
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          GO TO 100
        END IF
C                Control reaches this section if the error test has
C                failed MXFAIL or more times.  It is assumed that the
C                derivatives that have accumulated in the YH array have
C                errors of the wrong order.  Hence the first derivative
C                is recomputed, the order is set to 1, and the step is
C                retried.
        NFAIL = 0
        JTASK = 2
        DO 215 I = 1,N
 215      Y(I) = YH(I,1)
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              UROUND, USERS, Y, YWT,  H, MNTOLD, MTROLD, NFE, RC,
     8              YH,  A, CONVRG, EL, FAC, IER, IPVT, NQ, NWAIT, RH,
     8              RMAX, SAVE2, TQ, TREND, ISWFLG, JSTATE)
        RMAX = RMNORM
        IF (N .EQ. 0) GO TO 440
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (DBLE(NSTEP-1)*AVGH + H)/DBLE(NSTEP)
      AVGORD = (DBLE(NSTEP-1)*AVGORD + DBLE(NQ))/DBLE(NSTEP)
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
C
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.D0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.D0/DBLE(NQ+1)))
            IF (HS .GT. 1.2D0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.D0
              RMAX = RMNORM
              TREND = 1.D0
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.D0/DBLE(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.D0/DBLE(NQ+1)))
            HN = MIN(HN, 1.D0/(2.D0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.D0
              CONVRG = .FALSE.
              CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          END IF
        END IF
      END IF
      IF (SWITCH) THEN
        MINT = 2
        MNTOLD = MINT
        MITER = MTRSV
        MTROLD = MITER
        MAXORD = MIN(MXRDSV, 5)
        NQ = MIN(NQ, MAXORD)
        RC = 0.D0
        RMAX = RMNORM
        TREND = 1.D0
        CALL DDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C                           Consider changing H if NWAIT = 1.  Otherwise
C                           decrease NWAIT by 1.  If NWAIT is then 1 and
C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
C                           a possible order increase on the next step.
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (RH.GT.TRSHLD) CALL DDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
      ELSE IF (NWAIT .GT. 1) THEN
        NWAIT = NWAIT - 1
        IF (NWAIT .EQ. 1 .AND. NQ .LT. MAXORD) THEN
          DO 250 I = 1,NDE
 250        YH(I,MAXORD+1) = SAVE1(I)
        END IF
C             If a change in H is considered, an increase or decrease in
C             order by one is considered also.  A change in H is made
C             only if it is by a factor of at least TRSHLD.  Factors
C             RH1, RH2, and RH3 are computed, by which H could be
C             multiplied at order NQ - 1, order NQ, or order NQ + 1,
C             respectively.  The largest of these is determined and the
C             new order chosen accordingly.  If the order is to be
C             increased, we compute one additional scaled derivative.
C             If there is a change of order, reset NQ and the
C             coefficients.  In any case H is reset according to RH and
C             the YH array is rescaled.
      ELSE
        IF (NQ .EQ. 1) THEN
          RH1 = 0.D0
        ELSE
          DO 270 I = 1,NDE
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ERDN = DNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(DBLE(NDE)))
          RH1 = 1.D0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.D0/DBLE(NQ)))
        END IF
        RH2 = 1.D0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.D0/DBLE(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.D0
        ELSE
          DO 290 I = 1,NDE
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ERUP = DNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(DBLE(NDE)))
          RH3 = 1.D0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.D0/DBLE(NQ+2)))
        END IF
        IF (RH1 .GT. RH2 .AND. RH1 .GE. RH3) THEN
          RH = RH1
          IF (RH .LE. TRSHLD) GO TO 380
          NQ = NQ - 1
          RC = RC*EL(1,NQ)/EL(1,NQ+1)
        ELSE IF (RH2 .GE. RH1 .AND. RH2 .GE. RH3) THEN
          RH = RH2
          IF (RH .LE. TRSHLD) GO TO 380
        ELSE
          RH = RH3
          IF (RH .LE. TRSHLD) GO TO 380
          DO 360 I = 1,N
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/DBLE(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.D0) RH = MIN(RH, 1.D0/(2.D0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
        RMAX = RMNORM
 380    NWAIT = NQ + 2
      END IF
C               All returns are made through this section.  H is saved
C               in HOLD to allow the caller to change H on the next step
      JSTATE = 1
      HOLD = H
      RETURN
C
 400  JSTATE = 2
      HOLD = H
      DO 405 I = 1,N
 405    Y(I) = YH(I,1)
      RETURN
C
 410  JSTATE = 3
      HOLD = H
      RETURN
C
 420  JSTATE = 4
      HOLD = H
      RETURN
C
 430  T = TOLD
      CALL DDPSC (-1, NSV, NQ,  YH)
      DO 435 I = 1,NSV
 435    Y(I) = YH(I,1)
 440  HOLD = H
      RETURN
      END
      SUBROUTINE DDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DDZRO
C***REFER TO  DDRIV3
C     This is a special purpose version of ZEROIN, modified for use with
C     the DDRIV1 package.
C
C     Sandia Mathematical Program Library
C     Mathematical Computing Services Division 5422
C     Sandia Laboratories
C     P. O. Box 5800
C     Albuquerque, New Mexico  87115
C     Control Data 6600 Version 4.5, 1 November 1971
C
C     ABSTRACT
C        ZEROIN searches for a zero of a function F(N, T, Y, IROOT)
C        between the given values B and C until the width of the
C        interval (B, C) has collapsed to within a tolerance specified
C        by the stopping criterion, ABS(B - C) .LE. 2.*(RW*ABS(B) + AE).
C
C     Description of parameters
C        F     - Name of the external function, which returns a
C                double precision result.  This name must be in an
C                EXTERNAL statement in the calling program.
C        B     - One end of the interval (B, C).  The value returned for
C                B usually is the better approximation to a zero of F.
C        C     - The other end of the interval (B, C).
C        RE    - Relative error used for RW in the stopping criterion.
C                If the requested RE is less than machine precision,
C                then RW is set to approximately machine precision.
C        AE    - Absolute error used in the stopping criterion.  If the
C                given interval (B, C) contains the origin, then a
C                nonzero value should be chosen for AE.
C
C     REFERENCES
C       1.  L F Shampine and H A Watts, ZEROIN, A Root-Solving Routine,
C           SC-TM-70-631, Sept 1970.
C       2.  T J Dekker, Finding a Zero by Means of Successive Linear
C           Interpolation, "Constructive Aspects of the Fundamental
C           Theorem of Algebra", edited by B Dejon and P Henrici, 1969.
C***ROUTINES CALLED  DDNTP
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  870511   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDZRO
      DOUBLE PRECISION A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDZRO
      ER = 4.D0*UROUND
      RW = MAX(RE, ER)
      IC = 0
      ACBS = ABS(B - C)
      A = C
      FA = FC
      KOUNT = 0
C                                                    Perform interchange
 10   IF (ABS(FC) .LT. ABS(FB)) THEN
        A = B
        FA = FB
        B = C
        FB = FC
        C = A
        FC = FA
      END IF
      CMB = 0.5D0*(C - B)
      ACMB = ABS(CMB)
      TOL = RW*ABS(B) + AE
C                                                Test stopping criterion
      IF (ACMB .LE. TOL) RETURN
      IF (KOUNT .GT. 50) RETURN
C                                    Calculate new iterate implicitly as
C                                    B + P/Q, where we arrange P .GE. 0.
C                         The implicit form is used to prevent overflow.
      P = (B - A)*FB
      Q = FA - FB
      IF (P .LT. 0.D0) THEN
        P = -P
        Q = -Q
      END IF
C                          Update A and check for satisfactory reduction
C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.D0*ACMB .GE. ACBS) THEN
C                                                                 Bisect
          B = 0.5D0*(C + B)
          GO TO 20
        END IF
        IC = 0
      END IF
      ACBS = ACMB
C                                            Test for too small a change
      IF (P .LE. ABS(Q)*TOL) THEN
C                                                 Increment by tolerance
        B = B + SIGN(TOL, CMB)
C                                               Root ought to be between
C                                               B and (C + B)/2.
      ELSE IF (P .LT. CMB*Q) THEN
C                                                            Interpolate
        B = B + P/Q
      ELSE
C                                                                 Bisect
        B = 0.5D0*(C + B)
      END IF
C                                             Have completed computation
C                                             for new iterate B.
 20   CALL DDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (N .EQ. 0) RETURN
      IF (FB .EQ. 0.D0) RETURN
      KOUNT = KOUNT + 1
C
C             Decide whether next step is interpolation or extrapolation
C
      IF (SIGN(1.0D0, FB) .EQ. SIGN(1.0D0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
      END
      SUBROUTINE DEA(NEWFLG,SVALUE,LIMEXP,RESULT,ABSERR,EPSTAB,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DEA
C***DATE WRITTEN   800101  (YYMMDD)
C***REVISION DATE  880908  (YYMMDD)
C***CATEGORY NO.  E5
C***KEYWORDS  CONVERGENCE ACCELERATION,EPSILON ALGORITHM,EXTRAPOLATION
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONKER-KAPENGA, ELISE,WESTERN MICHIGAN UNIVERSITY
C           KAHANER, DAVID K., NATIONAL BUREAU OF STANDARDS
C           STARKENBURG, C. B., NATIONAL BUREAU OF STANDARDS
C***PURPOSE  Given a slowly convergent sequence, this routine attempts
C            to extrapolate nonlinearly to a better estimate of the
C            sequence's limiting value, thus improving the rate of
C            convergence. Routine is based on the epsilon algorithm
C            of P. Wynn. An estimate of the absolute error is also
C            given.
C***DESCRIPTION
C
C              Epsilon algorithm. Standard fortran subroutine.
C              Double precision version.
C
C       A R G U M E N T S   I N   T H E   C A L L   S E Q U E N C E
C
C              NEWFLG - LOGICAL                       (INPUT and OUTPUT)
C                       On the first call to DEA set NEWFLG to .TRUE.
C                       (indicating a new sequence). DEA will set NEWFLG
C                       to .FALSE.
C
C              SVALUE - DOUBLE PRECISION                         (INPUT)
C                       On the first call to DEA set SVALUE to the first
C                       term in the sequence. On subsequent calls set
C                       SVALUE to the subsequent sequence value.
C
C              LIMEXP - INTEGER                                  (INPUT)
C                       An integer equal to or greater than the total
C                       number of sequence terms to be evaluated. Do not
C                       change the value of LIMEXP until a new sequence
C                       is evaluated (NEWFLG=.TRUE.).  LIMEXP .GE. 3
C
C              RESULT - DOUBLE PRECISION                        (OUTPUT)
C                       Best approximation to the sequence's limit.
C
C              ABSERR - DOUBLE PRECISION                        (OUTPUT)
C                       Estimate of the absolute error.
C
C              EPSTAB - DOUBLE PRECISION                        (OUTPUT)
C                       Workvector of DIMENSION at least (LIMEXP+7).
C
C              IERR   - INTEGER                                 (OUTPUT)
C                       IERR=0 Normal termination of the routine.
C                       IERR=1 The input is invalid because LIMEXP.LT.3.
C
C    T Y P I C A L   P R O B L E M   S E T U P
C
C   This sample problem uses the trapezoidal rule to evaluate the
C   integral of the sin function from 0.0 to 0.5*PI (value = 1.0). The
C   program implements the trapezoidal rule 8 times creating an
C   increasingly accurate sequence of approximations to the integral.
C   Each time the trapezoidal rule is used, it uses twice as many
C   panels as the time before. DEA is called to obtain even more
C   accurate estimates.
C
C      PROGRAM SAMPLE
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      DOUBLE PRECISION EPSTAB(57)
CC                                     [57 = LIMEXP + 7]
C      LOGICAL NEWFLG
C      EXTERNAL F
C      DATA LIMEXP/50/
C      WRITE(*,*) ' NO. PANELS          TRAP. APPROX'
C     *           ,'        APPROX W/EA         ABSERR'
C      WRITE(*,*)
C      HALFPI = DASIN(1.0D+00)
CC                                     [UPPER INTEGRATION LIMIT = PI/2]
C      NEWFLG = .TRUE.
CC                                     [SET FLAG - 1ST DEA CALL]
C      DO 10 I = 0,7
C         NPARTS = 2 ** I
C         WIDTH = HALFPI/NPARTS
C         APPROX = 0.5D+00 * WIDTH * (F(0.0D+00) + F(HALFPI))
C         DO 11 J = 1,NPARTS-1
C            APPROX = APPROX + F(J * WIDTH) * WIDTH
C  11     CONTINUE
CC                                     [END TRAPEZOIDAL RULE APPROX]
C         SVALUE = APPROX
CC                                     [SVALUE = NEW SEQUENCE VALUE]
C         CALL DEA(NEWFLG,SVALUE,LIMEXP,RESULT,ABSERR,EPSTAB,IERR)
CC                                     [CALL DEA FOR BETTER ESTIMATE]
C         WRITE(*,12) NPARTS,APPROX,RESULT,ABSERR
C  12     FORMAT('   ',I4,T20,F16.13,T40,F16.13,T60,D11.4)
C  10  CONTINUE
C      STOP
C      END
C
C      DOUBLE PRECISION FUNCTION F(X)
C      DOUBLE PRECISION X
C      F = DSIN(X)
CC                                     [INTEGRAND]
C      RETURN
C      END
C
C   Output from the above program will be:
C
C  NO. PANELS          TRAP. APPROX        APPROX W/EA         ABSERR
C
C      1              .7853981633974      .7853981633974      .7854D+00
C      2              .9480594489685      .9480594489685      .9760D+00
C      4              .9871158009728      .9994567212570      .2141D+00
C      8              .9967851718862      .9999667417647      .5100D-03
C     16              .9991966804851      .9999998781041      .5763D-03
C     32              .9997991943200      .9999999981026      .5767D-03
C     64              .9999498000921      .9999999999982      .3338D-04
C    128              .9999874501175     1.0000000000000      .1238D-06
C
C-----------------------------------------------------------------------
C***REFERENCES  "Acceleration de la convergence en analyse numerique",
C                 C. Brezinski, "Lecture Notes in Math.", vol. 584,
C                 Springer-Verlag, New York, 1977.
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  DEA
      DOUBLE PRECISION ABSERR,DELTA1,DELTA2,DELTA3,DRELPR,DEPRN,
     1   D1MACH,EPSTAB(*),ERROR,ERR1,ERR2,ERR3,E0,E1,E2,E3,RES,
     2   RESULT,RES3LA(3),SS,SVALUE,TOL1,TOL2,TOL3
      INTEGER I,IB,IB2,IE,IERR,IN,K1,K2,K3,LIMEXP,N,NEWELM,NUM,NRES
      LOGICAL NEWFLG
C
C
C           LIMEXP  is the maximum number of elements the
C           epsilon table data can contain. The epsilon table
C           is stored in the first (LIMEXP+2) entries of EPSTAB.
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           E0,E1,E2,E3 - DOUBLE PRECISION
C                         The 4 elements on which the computation of
C                         a new element in the epsilon table is based.
C           NRES   - INTEGER
C                    Number of extrapolation results actually
C                    generated by the epsilon algorithm in prior
C                    calls to the routine.
C           NEWELM - INTEGER
C                    Number of elements to be computed in the
C                    new diagonal of the epsilon table. The
C                    condensed epsilon table is computed. Only
C                    those elements needed for the computation of
C                    the next diagonal are preserved.
C           RES    - DOUBLE PRECISION
C                    New element in the new diagonal of the
C                    epsilon table.
C           ERROR  - DOUBLE PRECISION
C                    An estimate of the absolute error of RES.
C                    Routine decides whether RESULT=RES or
C                    RESULT=SVALUE by comparing ERROR with
C                    ABSERR from the previous call.
C           RES3LA - DOUBLE PRECISION
C                    Vector of DIMENSION 3 containing at most
C                    the last 3 results.
C
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C            DRELPR  is the largest relative spacing.
C
C***FIRST EXECUTABLE STATEMENT  DEA
      IF(LIMEXP.LT.3) THEN
        IERR = 1
        CALL XERROR('LIMEXP IS LESS THAN 3',21,1,1)
        GO TO 110
      ENDIF
      IERR = 0
      RES3LA(1)=EPSTAB(LIMEXP+5)
      RES3LA(2)=EPSTAB(LIMEXP+6)
      RES3LA(3)=EPSTAB(LIMEXP+7)
      RESULT=SVALUE
      IF(NEWFLG) THEN
        N=1
        NRES=0
        NEWFLG=.FALSE.
        EPSTAB(N)=SVALUE
        ABSERR=ABS(RESULT)
        GO TO 100
      ELSE
        N = INT(EPSTAB(LIMEXP+3))
        NRES=INT(EPSTAB(LIMEXP+4))
        IF(N.EQ.2) THEN
          EPSTAB(N)=SVALUE
          ABSERR=.6D+01*ABS(RESULT-EPSTAB(1))
          GO TO 100
        ENDIF
      ENDIF
      EPSTAB(N)=SVALUE
      DRELPR=D1MACH(4)
      DEPRN=1.0D+01*DRELPR
      EPSTAB(N+2)=EPSTAB(N)
      NEWELM=(N-1)/2
      NUM=N
      K1=N
      DO 40 I=1,NEWELM
        K2=K1-1
        K3=K1-2
        RES=EPSTAB(K1+2)
        E0=EPSTAB(K3)
        E1=EPSTAB(K2)
        E2=RES
        DELTA2=E2-E1
        ERR2=ABS(DELTA2)
        TOL2=MAX(ABS(E2),ABS(E1))*DRELPR
        DELTA3=E1-E0
        ERR3=ABS(DELTA3)
        TOL3=MAX(ABS(E1),ABS(E0))*DRELPR
        IF(ERR2.GT.TOL2.OR.ERR3.GT.TOL3) GO TO 10
C
C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
C           ACCURACY, CONVERGENCE IS ASSUMED.
C           RESULT=E2
C           ABSERR=ABS(E1-E0)+ABS(E2-E1)
C
        RESULT=RES
        ABSERR=ERR2+ERR3
        GO TO 50
   10   IF(I.NE.1) THEN
          E3=EPSTAB(K1)
          EPSTAB(K1)=E1
          DELTA1=E1-E3
          ERR1=ABS(DELTA1)
          TOL1=MAX(ABS(E1),ABS(E3))*DRELPR
C
C           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
C           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
C
          IF(ERR1.LE.TOL1.OR.ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
          SS=0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
        ELSE
          EPSTAB(K1)=E1
          IF(ERR2.LE.TOL2.OR.ERR3.LE.TOL3) GO TO 20
          SS=0.1D+01/DELTA2-0.1D+01/DELTA3
        ENDIF
C
C           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
C           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
C           OF N
C
        IF(ABS(SS*E1).GT.0.1D-03) GO TO 30
   20   N=I+I-1
        IF(NRES.EQ.0) THEN
          ABSERR=ERR2+ERR3
          RESULT=RES
        ELSE IF(NRES.EQ.1) THEN
          RESULT=RES3LA(1)
        ELSE IF(NRES.EQ.2) THEN
          RESULT=RES3LA(2)
        ELSE
          RESULT=RES3LA(3)
        ENDIF
        GO TO 50
C
C           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
C           THE VALUE OF RESULT
C
   30   RES=E1+0.1D+01/SS
        EPSTAB(K1)=RES
        K1=K1-2
        ERROR=ERR2+ABS(RES-E2)+ERR3
        IF(NRES.EQ.0) GO TO 35
        IF(ERROR.GT.0.1E3*ABSERR) GO TO 40
   35   ABSERR=ERROR
        RESULT=RES
   40   CONTINUE
C
C           COMPUTE ERROR ESTIMATE
C
        IF(NRES.EQ.1) THEN
          ABSERR=ABS(RESULT-RES3LA(1))
        ELSE IF(NRES.EQ.2) THEN
          ABSERR=ABS(RESULT-RES3LA(2))+ABS(RESULT-RES3LA(1))
        ELSE IF(NRES.GT.2) THEN
          ABSERR=ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2))
     1          +ABS(RESULT-RES3LA(1))
        ENDIF
C
C           SHIFT THE TABLE
C
   50 IF(N.EQ.LIMEXP) N=2*(LIMEXP/2)-1
      IB=1
      IF((NUM/2)*2.EQ.NUM) IB=2
      IE=NEWELM+1
      DO 60 I=1,IE
        IB2=IB+2
        EPSTAB(IB)=EPSTAB(IB2)
        IB=IB2
   60 CONTINUE
      IF(NUM.EQ.N) GO TO 80
      IN=NUM-N+1
      DO 70 I=1,N
        EPSTAB(I)=EPSTAB(IN)
        IN=IN+1
   70 CONTINUE
C
C           UPDATE RES3LA
C
   80 IF(NRES.EQ.0) THEN
        RES3LA(1)=RESULT
      ELSE IF(NRES.EQ.1) THEN
        RES3LA(2)=RESULT
      ELSE IF(NRES.EQ.2) THEN
        RES3LA(3)=RESULT
      ELSE
        RES3LA(1)=RES3LA(2)
        RES3LA(2)=RES3LA(3)
        RES3LA(3)=RESULT
      ENDIF
   90 ABSERR=MAX(ABSERR,DEPRN*ABS(RESULT))
      NRES=NRES+1
  100 N=N+1
      EPSTAB(LIMEXP+3)=DBLE(N)
      EPSTAB(LIMEXP+4)=DBLE(NRES)
      EPSTAB(LIMEXP+5)=RES3LA(1)
      EPSTAB(LIMEXP+6)=RES3LA(2)
      EPSTAB(LIMEXP+7)=RES3LA(3)
  110 RETURN
      END
      DOUBLE PRECISION FUNCTION DENORM(N,X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DENORM
C***REFER TO  DNSQ,DNSQE
C
C        ******
C     FUNCTION DENORM
C
C     Given an N-vector X, this function calculates the
C     Euclidean norm of X.
C
C     The Euclidean norm is computed by accumulating the sum of
C     squares in three different sums. The sums of squares for the
C     small and large components are scaled so that no overflows
C     occur. Non-destructive underflows are permitted. Underflows
C     and overflows do not occur in the computation of the unscaled
C     sum of squares for the intermediate components.
C     The definitions of small, intermediate and large components
C     depend on two constants, RDWARF and RGIANT. The main
C     restrictions on these constants are that RDWARF**2 not
C     underflow and RGIANT**2 not overflow. The constants
C     given here are suitable for every known computer.
C
C     The function statement is
C
C       DOUBLE PRECISION FUNCTION DENORM(N,X)
C
C     where
C
C       N is a positive integer input variable.
C
C       X is an input array of length N.
C
C     Subprograms called
C
C       FORTRAN-Supplied ... DABS,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DENORM
      DOUBLE PRECISION DABS, DSQRT
      INTEGER I, N
      DOUBLE PRECISION AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3,
     1     X(N), X1MAX, X3MAX, XABS, ZERO
      SAVE ONE, ZERO, RDWARF, RGIANT
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
C***FIRST EXECUTABLE STATEMENT  DENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = DABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
C
C              SUM FOR LARGE COMPONENTS.
C
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
C
C              SUM FOR SMALL COMPONENTS.
C
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
C
C           SUM FOR INTERMEDIATE COMPONENTS.
C
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
C
C     CALCULATION OF NORM.
C
      IF (S1 .EQ. ZERO) GO TO 100
         DENORM = X1MAX*DSQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX)
     1         DENORM = DSQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX)
     1         DENORM = DSQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            DENORM = X3MAX*DSQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DERF(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DERF
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  840425   (YYMMDD)
C***CATEGORY NO.  C8A,L5A1E
C***KEYWORDS  DOUBLE PRECISION,ERROR FUNCTION,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. error function, ERF, of X.
C***DESCRIPTION
C
C DERF(X) calculates the double precision error function for double
C precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted error   1.28E-32
C                                         log weighted error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,DERFC,INITDS
C***END PROLOGUE  DERF
      DOUBLE PRECISION X, ERFCS(21), SQEPS, SQRTPI, XBIG, Y, D1MACH,
     1  DCSEVL, DERFC
      DATA ERF CS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERF CS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERF CS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERF CS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERF CS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERF CS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERF CS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERF CS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERF CS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERF CS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERF CS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERF CS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERF CS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERF CS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERF CS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERF CS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERF CS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERF CS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERF CS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERF CS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERF CS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA NTERF, XBIG, SQEPS / 0, 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DERF
      IF (NTERF.NE.0) GO TO 10
      NTERF = INITDS (ERFCS, 21, 0.1*SNGL(D1MACH(3)))
      XBIG = DSQRT (-DLOG(SQRTPI*D1MACH(3)))
      SQEPS = DSQRT (2.0D0*D1MACH(3))
C
 10   Y = DABS(X)
      IF (Y.GT.1.D0) GO TO 20
C
C ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
C
      IF (Y.LE.SQEPS) DERF = 2.0D0*X/SQRTPI
      IF (Y.GT.SQEPS) DERF = X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
C
 20   IF (Y.LE.XBIG) DERF = DSIGN (1.0D0-DERFC(Y), X)
      IF (Y.GT.XBIG) DERF = DSIGN (1.0D0, X)
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DERFC(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DERFC
C***DATE WRITTEN   770701   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C8A,L5A1E
C***KEYWORDS  COMPLEMENTARY ERROR FUNCTION,DOUBLE PRECISION,
C             SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. complementary error function, ERFC.
C***DESCRIPTION
C
C DERFC(X) calculates the double precision complementary error function
C for double precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted Error   1.28E-32
C                                         log weighted Error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00
C                                        with weighted Error   2.67E-32
C                                         log weighted Error  31.57
C                               significant figures required  30.31
C                                    decimal places required  32.42
C
C Series for ERFC       on the interval  0.          to  2.50000E-01
C                                        with weighted error   1.53E-31
C                                         log weighted error  30.82
C                               significant figures required  29.47
C                                    decimal places required  31.70
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,DCSEVL,INITDS,XERROR
C***END PROLOGUE  DERFC
      DOUBLE PRECISION X, ERFCS(21), ERFCCS(59), ERC2CS(49), SQEPS,
     1  SQRTPI, XMAX, XSML, Y,  D1MACH, DCSEVL
      DATA ERF CS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERF CS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERF CS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERF CS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERF CS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERF CS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERF CS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERF CS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERF CS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERF CS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERF CS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERF CS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERF CS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERF CS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERF CS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERF CS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERF CS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERF CS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERF CS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERF CS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERF CS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA ERC2CS(  1) / -.6960134660 2309501127 3915082619 7 D-1      /
      DATA ERC2CS(  2) / -.4110133936 2620893489 8221208466 6 D-1      /
      DATA ERC2CS(  3) / +.3914495866 6896268815 6114370524 4 D-2      /
      DATA ERC2CS(  4) / -.4906395650 5489791612 8093545077 4 D-3      /
      DATA ERC2CS(  5) / +.7157479001 3770363807 6089414182 5 D-4      /
      DATA ERC2CS(  6) / -.1153071634 1312328338 0823284791 2 D-4      /
      DATA ERC2CS(  7) / +.1994670590 2019976350 5231486770 9 D-5      /
      DATA ERC2CS(  8) / -.3642666471 5992228739 3611843071 1 D-6      /
      DATA ERC2CS(  9) / +.6944372610 0050125899 3127721463 3 D-7      /
      DATA ERC2CS( 10) / -.1371220902 1043660195 3460514121 0 D-7      /
      DATA ERC2CS( 11) / +.2788389661 0071371319 6386034808 7 D-8      /
      DATA ERC2CS( 12) / -.5814164724 3311615518 6479105031 6 D-9      /
      DATA ERC2CS( 13) / +.1238920491 7527531811 8016881795 0 D-9      /
      DATA ERC2CS( 14) / -.2690639145 3067434323 9042493788 9 D-10     /
      DATA ERC2CS( 15) / +.5942614350 8479109824 4470968384 0 D-11     /
      DATA ERC2CS( 16) / -.1332386735 7581195792 8775442057 0 D-11     /
      DATA ERC2CS( 17) / +.3028046806 1771320171 7369724330 4 D-12     /
      DATA ERC2CS( 18) / -.6966648814 9410325887 9586758895 4 D-13     /
      DATA ERC2CS( 19) / +.1620854541 0539229698 1289322762 8 D-13     /
      DATA ERC2CS( 20) / -.3809934465 2504919998 7691305772 9 D-14     /
      DATA ERC2CS( 21) / +.9040487815 9788311493 6897101297 5 D-15     /
      DATA ERC2CS( 22) / -.2164006195 0896073478 0981204700 3 D-15     /
      DATA ERC2CS( 23) / +.5222102233 9958549846 0798024417 2 D-16     /
      DATA ERC2CS( 24) / -.1269729602 3645553363 7241552778 0 D-16     /
      DATA ERC2CS( 25) / +.3109145504 2761975838 3622741295 1 D-17     /
      DATA ERC2CS( 26) / -.7663762920 3203855240 0956671481 1 D-18     /
      DATA ERC2CS( 27) / +.1900819251 3627452025 3692973329 0 D-18     /
      DATA ERC2CS( 28) / -.4742207279 0690395452 2565599996 5 D-19     /
      DATA ERC2CS( 29) / +.1189649200 0765283828 8068307845 1 D-19     /
      DATA ERC2CS( 30) / -.3000035590 3257802568 4527131306 6 D-20     /
      DATA ERC2CS( 31) / +.7602993453 0432461730 1938527709 8 D-21     /
      DATA ERC2CS( 32) / -.1935909447 6068728815 6981104913 0 D-21     /
      DATA ERC2CS( 33) / +.4951399124 7733378810 0004238677 3 D-22     /
      DATA ERC2CS( 34) / -.1271807481 3363718796 0862198988 8 D-22     /
      DATA ERC2CS( 35) / +.3280049600 4695130433 1584165205 3 D-23     /
      DATA ERC2CS( 36) / -.8492320176 8228965689 2479242239 9 D-24     /
      DATA ERC2CS( 37) / +.2206917892 8075602235 1987998719 9 D-24     /
      DATA ERC2CS( 38) / -.5755617245 6965284983 1281950719 9 D-25     /
      DATA ERC2CS( 39) / +.1506191533 6392342503 5414405119 9 D-25     /
      DATA ERC2CS( 40) / -.3954502959 0187969531 0428569599 9 D-26     /
      DATA ERC2CS( 41) / +.1041529704 1515009799 8464505173 3 D-26     /
      DATA ERC2CS( 42) / -.2751487795 2787650794 5017890133 3 D-27     /
      DATA ERC2CS( 43) / +.7290058205 4975574089 9770368000 0 D-28     /
      DATA ERC2CS( 44) / -.1936939645 9159478040 7750109866 6 D-28     /
      DATA ERC2CS( 45) / +.5160357112 0514872983 7005482666 6 D-29     /
      DATA ERC2CS( 46) / -.1378419322 1930940993 8964480000 0 D-29     /
      DATA ERC2CS( 47) / +.3691326793 1070690422 5109333333 3 D-30     /
      DATA ERC2CS( 48) / -.9909389590 6243654206 5322666666 6 D-31     /
      DATA ERC2CS( 49) / +.2666491705 1953884133 2394666666 6 D-31     /
      DATA ERFCCS(  1) / +.7151793102 0292477450 3697709496 D-1        /
      DATA ERFCCS(  2) / -.2653243433 7606715755 8893386681 D-1        /
      DATA ERFCCS(  3) / +.1711153977 9208558833 2699194606 D-2        /
      DATA ERFCCS(  4) / -.1637516634 5851788416 3746404749 D-3        /
      DATA ERFCCS(  5) / +.1987129350 0552036499 5974806758 D-4        /
      DATA ERFCCS(  6) / -.2843712412 7665550875 0175183152 D-5        /
      DATA ERFCCS(  7) / +.4606161308 9631303696 9379968464 D-6        /
      DATA ERFCCS(  8) / -.8227753025 8792084205 7766536366 D-7        /
      DATA ERFCCS(  9) / +.1592141872 7709011298 9358340826 D-7        /
      DATA ERFCCS( 10) / -.3295071362 2528432148 6631665072 D-8        /
      DATA ERFCCS( 11) / +.7223439760 4005554658 1261153890 D-9        /
      DATA ERFCCS( 12) / -.1664855813 3987295934 4695966886 D-9        /
      DATA ERFCCS( 13) / +.4010392588 2376648207 7671768814 D-10       /
      DATA ERFCCS( 14) / -.1004816214 4257311327 2170176283 D-10       /
      DATA ERFCCS( 15) / +.2608275913 3003338085 9341009439 D-11       /
      DATA ERFCCS( 16) / -.6991110560 4040248655 7697812476 D-12       /
      DATA ERFCCS( 17) / +.1929492333 2617070862 4205749803 D-12       /
      DATA ERFCCS( 18) / -.5470131188 7543310649 0125085271 D-13       /
      DATA ERFCCS( 19) / +.1589663309 7626974483 9084032762 D-13       /
      DATA ERFCCS( 20) / -.4726893980 1975548392 0369584290 D-14       /
      DATA ERFCCS( 21) / +.1435873376 7849847867 2873997840 D-14       /
      DATA ERFCCS( 22) / -.4449510561 8173583941 7250062829 D-15       /
      DATA ERFCCS( 23) / +.1404810884 7682334373 7305537466 D-15       /
      DATA ERFCCS( 24) / -.4513818387 7642108962 5963281623 D-16       /
      DATA ERFCCS( 25) / +.1474521541 0451330778 7018713262 D-16       /
      DATA ERFCCS( 26) / -.4892621406 9457761543 6841552532 D-17       /
      DATA ERFCCS( 27) / +.1647612141 4106467389 5301522827 D-17       /
      DATA ERFCCS( 28) / -.5626817176 3294080929 9928521323 D-18       /
      DATA ERFCCS( 29) / +.1947443382 2320785142 9197867821 D-18       /
      DATA ERFCCS( 30) / -.6826305642 9484207295 6664144723 D-19       /
      DATA ERFCCS( 31) / +.2421988887 2986492401 8301125438 D-19       /
      DATA ERFCCS( 32) / -.8693414133 5030704256 3800861857 D-20       /
      DATA ERFCCS( 33) / +.3155180346 2280855712 2363401262 D-20       /
      DATA ERFCCS( 34) / -.1157372324 0496087426 1239486742 D-20       /
      DATA ERFCCS( 35) / +.4288947161 6056539462 3737097442 D-21       /
      DATA ERFCCS( 36) / -.1605030742 0576168500 5737770964 D-21       /
      DATA ERFCCS( 37) / +.6063298757 4538026449 5069923027 D-22       /
      DATA ERFCCS( 38) / -.2311404251 6979584909 8840801367 D-22       /
      DATA ERFCCS( 39) / +.8888778540 6618855255 4702955697 D-23       /
      DATA ERFCCS( 40) / -.3447260576 6513765223 0718495566 D-23       /
      DATA ERFCCS( 41) / +.1347865460 2069650682 7582774181 D-23       /
      DATA ERFCCS( 42) / -.5311794071 1250217364 5873201807 D-24       /
      DATA ERFCCS( 43) / +.2109341058 6197831682 8954734537 D-24       /
      DATA ERFCCS( 44) / -.8438365587 9237891159 8133256738 D-25       /
      DATA ERFCCS( 45) / +.3399982524 9452089062 7359576337 D-25       /
      DATA ERFCCS( 46) / -.1379452388 0732420900 2238377110 D-25       /
      DATA ERFCCS( 47) / +.5634490311 8332526151 3392634811 D-26       /
      DATA ERFCCS( 48) / -.2316490434 4770654482 3427752700 D-26       /
      DATA ERFCCS( 49) / +.9584462844 6018101526 3158381226 D-27       /
      DATA ERFCCS( 50) / -.3990722880 3301097262 4224850193 D-27       /
      DATA ERFCCS( 51) / +.1672129225 9444773601 7228709669 D-27       /
      DATA ERFCCS( 52) / -.7045991522 7660138563 8803782587 D-28       /
      DATA ERFCCS( 53) / +.2979768402 8642063541 2357989444 D-28       /
      DATA ERFCCS( 54) / -.1262522466 4606192972 2422632994 D-28       /
      DATA ERFCCS( 55) / +.5395438704 5424879398 5299653154 D-29       /
      DATA ERFCCS( 56) / -.2380992882 5314591867 5346190062 D-29       /
      DATA ERFCCS( 57) / +.1099052830 1027615735 9726683750 D-29       /
      DATA ERFCCS( 58) / -.4867713741 6449657273 2518677435 D-30       /
      DATA ERFCCS( 59) / +.1525877264 1103575676 3200828211 D-30       /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA NTERF, NTERFC, NTERC2, XSML, XMAX, SQEPS / 3*0, 3*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DERFC
      IF (NTERF.NE.0) GO TO 10
      ETA = 0.1*SNGL(D1MACH(3))
      NTERF = INITDS (ERFCS, 21, ETA)
      NTERFC = INITDS (ERFCCS, 59, ETA)
      NTERC2 = INITDS (ERC2CS, 49, ETA)
C
      XSML = -DSQRT (-DLOG(SQRTPI*D1MACH(3)))
      XMAX = DSQRT (-DLOG(SQRTPI*D1MACH(1)) )
      XMAX = XMAX - 0.5D0*DLOG(XMAX)/XMAX - 0.01D0
      SQEPS = DSQRT (2.0D0*D1MACH(3))
C
 10   IF (X.GT.XSML) GO TO 20
C
C ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
C
      DERFC = 2.0D0
      RETURN
C
 20   IF (X.GT.XMAX) GO TO 40
      Y = DABS(X)
      IF (Y.GT.1.0D0) GO TO 30
C
C ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
C
      IF (Y.LT.SQEPS) DERFC = 1.0D0 - 2.0D0*X/SQRTPI
      IF (Y.GE.SQEPS) DERFC = 1.0D0 - X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
C
 30   Y = Y*Y
      IF (Y.LE.4.D0) DERFC = DEXP(-Y)/DABS(X) * (0.5D0 + DCSEVL (
     1  (8.D0/Y-5.D0)/3.D0, ERC2CS, NTERC2) )
      IF (Y.GT.4.D0) DERFC = DEXP(-Y)/DABS(X) * (0.5D0 + DCSEVL (
     1  8.D0/Y-1.D0, ERFCCS, NTERFC) )
      IF (X.LT.0.D0) DERFC = 2.0D0 - DERFC
      RETURN
C
 40   CALL XERROR ( 'DERFC   X SO BIG ERFC UNDERFLOWS', 32, 1, 1)
      DERFC = 0.D0

      RETURN
      END
      SUBROUTINE DEZFT1(N,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DEZFT1
C***REFER TO  DEZFTI
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DEZFT1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  DEZFT1
      TPI = 8.0D0*ATAN(1.0D0)
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = TPI/DBLE(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 111 K1=1,NFM1
         IP = IFAC(K1+2)
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         ARG1 = DBLE(L1)*ARGH
         CH1 = 1.0D0
         SH1 = 0.0D0
         DCH1 = COS(ARG1)
         DSH1 = SIN(ARG1)
         DO 110 J=1,IPM
            CH1H = DCH1*CH1-DSH1*SH1
            SH1 = DCH1*SH1+DSH1*CH1
            CH1 = CH1H
            I = IS+2
            WA(I-1) = CH1
            WA(I) = SH1
            IF (IDO .LT. 5) GO TO 109
            DO 108 II=5,IDO,2
               I = I+2
               WA(I-1) = CH1*WA(I-3)-SH1*WA(I-2)
               WA(I) = CH1*WA(I-2)+SH1*WA(I-3)
  108       CONTINUE
  109       IS = IS+IDO
  110    CONTINUE
         L1 = L2
  111 CONTINUE
      RETURN
      END
      SUBROUTINE DEZFTB(N,R,AZERO,A,B,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DEZFTB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  A Simplified double precision, periodic,
C              backward transform
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C  Subroutine DEZFTB computes a d. p. periodic sequence from its
C  Fourier coefficients (Fourier synthesis).  The transform is
C  defined below at Output Parameter R.  DEZFTB is a simplified
C  but slower version of DRFFTB.
C
C  Input Parameters
C
C  N       the length of the output array R.  The method is most
C          efficient when N is the product of small primes.
C
C  AZERO   the constant Fourier coefficient
C
C  A,B     arrays which contain the remaining Fourier coefficients.
C          These arrays are not destroyed.
C
C          The length of these arrays depends on whether N is even or
C          odd.
C
C          If N is even, N/2    locations are required.
C          If N is odd, (N-1)/2 locations are required
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls DEZFTB.  The WSAVE array must be
C          initialized by calling subroutine DEZFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by EZFFTF and DEZFTB.
C
C
C  Output Parameters
C
C  R       if N is even, define KMAX=N/2
C          if N is odd,  define KMAX=(N-1)/2
C
C          Then for I=1,...,N
C
C               R(I)=AZERO plus the sum from K=1 to K=KMAX of
C
C               A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
C
C  ********************* Complex Notation **************************
C
C          For J=1,...,N
C
C          R(J) equals the sum from K=-KMAX to K=KMAX of
C
C               C(K)*EXP(I*K*(J-1)*2*PI/N)
C
C          where
C
C               C(K) = .5D0*CMPLX(A(K),-B(K))   for K=1,...,KMAX
C
C               C(-K) = CONJG(C(K))
C
C               C(0) = AZERO
C
C                    and I=SQRT(-1)
C
C  *************** Amplitude - Phase Notation ***********************
C
C          For I=1,...,N
C
C          R(I) equals AZERO plus the sum from K=1 to K=KMAX of
C
C               ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
C
C          where
C
C               ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
C
C               COS(BETA(K))=A(K)/ALPHA(K)
C
C               SIN(BETA(K))=-B(K)/ALPHA(K)
C
C  *                                                                   *
C  *   References                                                      *
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DRFFTB
C***END PROLOGUE  DEZFTB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       R(*)       ,A(*)       ,B(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DEZFTB
      IF (N-2) 101,102,103
  101 R(1) = AZERO
      RETURN
  102 R(1) = AZERO+A(1)
      R(2) = AZERO-A(1)
      RETURN
  103 NS2 = (N-1)/2
      DO 104 I=1,NS2
         R(2*I) = .5D0*A(I)
         R(2*I+1) = -.5D0*B(I)
  104 CONTINUE
      R(1) = AZERO
      IF (MOD(N,2) .EQ. 0) R(N) = A(NS2+1)
      CALL DRFFTB (N,R,WSAVE(N+1))
      RETURN
      END
      SUBROUTINE DEZFTF(N,R,AZERO,A,B,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DEZFTF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  A simplified d. p., periodic, forward transform
C***DESCRIPTION
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C  Subroutine DEZFTF computes the Fourier coefficients of a d. p.
C  perodic sequence (Fourier analysis).  The transform is defined
C  below at Output Parameters AZERO, A and B.  DEZFTF is a simplified
C  but slower version of DRFFTF.
C
C  Input Parameters
C
C  N       the length of the array R to be transformed.  The method
C          is must efficient when N is the product of small primes.
C
C  R       a d. p. array of length N which contains the sequence
C          to be transformed.  R is not destroyed.
C
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15
C          in the program that calls DEZFTF.  The WSAVE array must be
C          initialized by calling subroutine DEZFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by DEZFTF and EZFFTB.
C
C  Output Parameters
C
C  AZERO   the sum from I=1 to I=N of R(I)/N
C
C  A,B     for N even B(N/2)=0.0D0 and A(N/2) is the sum from I=1 to
C          I=N of (-1)**(I-1)*R(I)/N
C
C          for N even define KMAX=N/2-1
C          for N odd  define KMAX=(N-1)/2
C
C          then for  k=1,...,KMAX
C
C               A(K) equals the sum from I=1 to I=N of
C
C                    2.0D0/N*R(I)*COS(K*(I-1)*2*PI/N)
C
C               B(K) equals the sum from I=1 to I=N of
C
C                    2.0D0/N*R(I)*SIN(K*(I-1)*2*PI/N)
C
C  *                                                                   *
C  *   References                                                      *
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DRFFTF
C***END PROLOGUE  DEZFTF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       R(*)       ,A(*)       ,B(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DEZFTF
      IF (N-2) 101,102,103
  101 AZERO = R(1)
      RETURN
  102 AZERO = 0.5D0*(R(1)+R(2))
      A(1) = 0.5D0*(R(1)-R(2))
      RETURN
  103 DO 104 I=1,N
         WSAVE(I) = R(I)
  104 CONTINUE
      CALL DRFFTF (N,WSAVE,WSAVE(N+1))
      CF = 2.0D0/DBLE(N)
      CFM = -CF
      AZERO = 0.5D0*CF*WSAVE(1)
      NS2 = (N+1)/2
      NS2M = NS2-1
      DO 105 I=1,NS2M
         A(I) = CF*WSAVE(2*I)
         B(I) = CFM*WSAVE(2*I+1)
  105 CONTINUE
      IF (MOD(N,2) .EQ. 0) THEN
          A(NS2) = 0.5D0*CF*WSAVE(N)
          B(NS2) = 0.0D0
      ENDIF
      RETURN
      END
      SUBROUTINE DEZFTI(N,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DEZFTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize DEZFTF and DEZFTB
C***DESCRIPTION
C
C  Subroutine DEZFTI initializes the array WSAVE which is used in
C  both DEZFTF and DEZFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least 3*N+15.
C          The same work array can be used for both DEZFTF and DEZFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DEZFT1
C***END PROLOGUE  DEZFTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DEZFTI
      IF (N .EQ. 1) RETURN
      CALL DEZFT1 (N,WSAVE(2*N+1),WSAVE(3*N+1))
      RETURN
      END
      SUBROUTINE DFALTD(N,X,TYPSIZ,FSCALE,METHOD,IEXP,MSG,NDIGIT,
     +     ITNLIM,IAGFLG,IAHFLG,IPR,DLT,GRADTL,STEPMX,STEPTL)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C SET DEFAULT VALUES FOR EACH INPUT VARIABLE TO
C MINIMIZATION ALGORITHM.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> INITIAL GUESS TO SOLUTION (TO COMPUTE MAX STEP SIZE)
C TYPSIZ(N)   <--  TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE      <--  ESTIMATE OF SCALE OF MINIMIZATION FUNCTION
C METHOD      <--  ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C IEXP        <--  =0 IF MINIMIZATION FUNCTION NOT EXPENSIVE TO EVALUATE
C MSG         <--  MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT
C NDIGIT      <--  NUMBER OF GOOD DIGITS IN MINIMIZATION FUNCTION
C ITNLIM      <--  MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG      <--  =0 IF ANALYTIC GRADIENT NOT SUPPLIED
C IAHFLG      <--  =0 IF ANALYTIC HESSIAN NOT SUPPLIED
C IPR         <--  DEVICE TO WHICH TO SEND OUTPUT
C DLT         <--  TRUST REGION RADIUS
C GRADTL      <--  TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE ENOUGH
C                  TO ZERO TO TERMINATE ALGORITHM
C STEPMX      <--  VALUE OF ZERO TO TRIP DEFAULT MAXIMUM IN OPTCHD
C STEPTL      <--  TOLERANCE AT WHICH SUCCESSIVE ITERATES CONSIDERED
C                  CLOSE ENOUGH TO TERMINATE ALGORITHM
C
      DIMENSION TYPSIZ(N),X(N)
      X(N)=X(N)
C
C SET TYPICAL SIZE OF X AND MINIMIZATION FUNCTION
      DO 10 I=1,N
        TYPSIZ(I)=1.0D0
   10 CONTINUE
      FSCALE=1.0D0
C
C SET TOLERANCES
      DLT=-1.0D0
      EPSM=D1MACH(4)
      GRADTL=EPSM**(1.0D0/3.0D0)
      STEPMX=0.0D0
      STEPTL=SQRT(EPSM)
C
C SET FLAGS
      METHOD=1
      IEXP=1
C---
C     MSG=0
      MSG=9
      NDIGIT=-1
      ITNLIM=150
      IAGFLG=0
      IAHFLG=0
      IPR=I1MACH(2)

      RETURN
      END
      SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
     1   WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DFDJC1
C***REFER TO  DNSQ,DNSQE
C
C      *******
C     SUBROUTINE DFDJC1
C
C     This subroutine computes a forward-difference approximation
C     to the N by N Jacobian matrix associated with a specified
C     problem of N functions in N variables. If the Jacobian has
C     a banded form, then function evaluations are saved by only
C     approximating the nonzero terms.
C
C     The subroutine statement is
C
C       SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
C                         WA1,WA2)
C
C     where
C
C       FCN is the name of the user-supplied subroutine which
C         calculates the functions. FCN must be declared
C         in an EXTERNAL statement in the user calling
C         program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         Calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C
C         The value of IFLAG should not be changed by FCN unless
C         the user wants to terminate execution of DFDJC1.
C         In this case set IFLAG to a negative integer.
C
C       N is a positive integer input variable set to the number
C         of functions and variables.
C
C       X is an input array of length N.
C
C       FVEC is an input array of length N which must contain the
C         functions evaluated at X.
C
C       FJAC is an output N by N array which contains the
C         approximation to the Jacobian matrix evaluated at X.
C
C       LDFJAC is a positive integer input variable not less than N
C         which specifies the leading dimension of the array FJAC.
C
C       IFLAG is an integer variable which can be used to terminate
C         the execution of DFDJC1. See description of FCN.
C
C       ML is a nonnegative integer input variable which specifies
C         the number of subdiagonals within the band of the
C         Jacobian matrix. If the Jacobian is not banded, set
C         ML to at least N - 1.
C
C       EPSFCN is an input variable used in determining a suitable
C         step length for the forward-difference approximation. This
C         approximation assumes that the relative errors in the
C         functions are of the order of EPSFCN. If EPSFCN is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
C
C       MU is a nonnegative integer input variable which specifies
C         the number of superdiagonals within the band of the
C         Jacobian matrix. If the Jacobian is not banded, set
C         MU to at least N - 1.
C
C       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
C         least N, then the Jacobian is considered dense, and WA2 is
C         not referenced.
C
C     Subprograms called
C
C       SLATEC-supplied ... D1MACH
C
C       FORTRAN-supplied ... DABS,DMAX1,DSQRT
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DFDJC1
      DOUBLE PRECISION D1MACH
      DOUBLE PRECISION DABS, DMAX1, DSQRT
      INTEGER I, IFLAG, J, K, LDFJAC, ML, MSUM, MU, N
      DOUBLE PRECISION EPS, EPSFCN, EPSMCH, FJAC(LDFJAC,N),
     1     FVEC(N), H, TEMP, WA1(N), WA2(N), X(N), ZERO
      SAVE ZERO
      DATA ZERO /0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
C***FIRST EXECUTABLE STATEMENT  DFDJC1
      EPSMCH = D1MACH(4)
C
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
      MSUM = ML + MU + 1
      IF (MSUM .LT. N) GO TO 40
C
C        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
C
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*DABS(TEMP)
            IF (H .EQ. ZERO) H = EPS
            X(J) = TEMP + H
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
C
C        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
C
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*DABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*DABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I .GE. J - MU .AND. I .LE. J + ML)
     1               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
      DOUBLE PRECISION FUNCTION DFMIN(AX,BX,F,TOL)

c*********************************************************************72
c
      DOUBLE PRECISION AX,BX,F,TOL
C***BEGIN PROLOGUE  DFMIN
C***CATEGORY NO.    G1A2
C***KEYWORD(S)  ONE-DIMENSIONAL MINIMIZATION, UNIMODAL FUNCTION
C***AUTHOR  R. BRENT
C***DATE WRITTEN    730101  (YYMMDD)
C***PURPOSE
C     An approximation to the point where  F  attains a minimum  on
C     the interval (AX,BX) is determined as the value of the function
C     DFMIN.
C
C PARAMETERS
C
C INPUT
C
C  AX    (double precision)  left endpoint of initial interval
C  BX    (double precision) right endpoint of initial interval
C  F     function subprogram which evaluates F(X)  for any  X
C        in the interval  (AX,BX)
C  TOL   (double precision) desired length of the interval of uncertainty 
C        of the final result ( .ge. 0.0)
C
C
C OUTPUT
C
C DFMIN   abcissa approximating the minimizer of F
C AX     lower bound for minimizer
C BX     upper bound for minimizer
C
C***DESCRIPTION
C
C     The method used is a combination of golden section search and
C     successive parabolic interpolation.  Convergence is never much 
C     slower than that for a Fibonacci search.  If F has a continuous 
C     second derivative which is positive at the minimum (which is not
C     at AX or BX), then convergence is superlinear, and usually of the 
C     order of about 1.324....
C
C     The function F is never evaluated at two points closer together
C     than EPS*ABS(DFMIN) + (TOL/3), where EPS is approximately the 
C     square root of the relative machine precision.  If F is a unimodal
C     function and the computed values of F are always unimodal when
C     separated by at least EPS*ABS(XSTAR) + (TOL/3), then DFMIN 
C     approximates the abcissa of the global minimum of F on the 
C     interval AX,BX with an error less than 3*EPS*ABS(DFMIN) + TOL.  
C     If F is not unimodal, then DFMIN may approximate a local, but 
C     perhaps non-global, minimum to the same accuracy.
C
C     This function subprogram is a slightly modified version of the
C     ALGOL 60 procedure LOCALMIN given in Richard Brent, Algorithms for
C     Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
C
C***REFERENCE(S)
C     Richard Brent, Algorithms for Minimization Without Derivatives,
C     Prentice-Hall, Inc. (1973).
C***ROUTINES CALLED   NONE
C***END PROLOGUE
      DOUBLE PRECISION  A,B,C,D,E,EPS,XM,P,Q,R,TOL1,TOL2,U,V,W
      DOUBLE PRECISION  FU,FV,FW,FX,X
      DOUBLE PRECISION  DABS,DSQRT,DSIGN
C
C  C is the squared inverse of the golden ratio
C
      C = 0.5D0*(3. - DSQRT(5.0D0))
C
C  EPS is approximately the square root of the relative machine
C  precision.
C
      EPS = 1.0D00
   10 EPS = EPS/2.0D00
      TOL1 = 1.0D0 + EPS
      IF (TOL1 .GT. 1.0D00) GO TO 10
      EPS = DSQRT(EPS)
C
C  initialization
C
      A = AX
      B = BX
      V = A + C*(B - A)
      W = V
      X = V
      E = 0.0D0
      FX = F(X)
      FV = FX
      FW = FX
C
C  main loop starts here
C
   20 XM = 0.5D0*(A + B)
      TOL1 = EPS*DABS(X) + TOL/3.0D0
      TOL2 = 2.0D0*TOL1
C
C  check stopping criterion
C
      IF (DABS(X - XM) .LE. (TOL2 - 0.5D0*(B - A))) GO TO 90
C
C is golden-section necessary
C
      IF (DABS(E) .LE. TOL1) GO TO 40
C
C  fit parabola
C
      R = (X - W)*(FX - FV)
      Q = (X - V)*(FX - FW)
      P = (X - V)*Q - (X - W)*R
      Q = 2.0D00*(Q - R)
      IF (Q .GT. 0.0D0) P = -P
      Q =  DABS(Q)
      R = E
      E = D
C
C  is parabola acceptable
C
   30 IF (DABS(P) .GE. DABS(0.5D0*Q*R)) GO TO 40
      IF (P .LE. Q*(A - X)) GO TO 40
      IF (P .GE. Q*(B - X)) GO TO 40
C
C  a parabolic interpolation step
C
      D = P/Q
      U = X + D
C
C  F must not be evaluated too close to AX or BX
C
      IF ((U - A) .LT. TOL2) D = DSIGN(TOL1, XM - X)
      IF ((B - U) .LT. TOL2) D = DSIGN(TOL1, XM - X)
      GO TO 50
C
C  a golden-section step
C
   40 IF (X .GE. XM) E = A - X
      IF (X .LT. XM) E = B - X
      D = C*E
C
C  F must not be evaluated too close to X
C
   50 IF (DABS(D) .GE. TOL1) U = X + D
      IF (DABS(D) .LT. TOL1) U = X + DSIGN(TOL1, D)
      FU = F(U)
C
C  update  A, B, V, W, and X
C
      IF (FU .GT. FX) GO TO 60
      IF (U .GE. X) A = X
      IF (U .LT. X) B = X
      V = W
      FV = FW
      W = X
      FW = FX
      X = U
      FX = FU
      GO TO 20
   60 IF (U .LT. X) A = U
      IF (U .GE. X) B = U
      IF (FU .LE. FW) GO TO 70
      IF (W .EQ. X) GO TO 70
      IF (FU .LE. FV) GO TO 80
      IF (V .EQ. X) GO TO 80
      IF (V .EQ. W) GO TO 80
      GO TO 20
   70 V = W
      FV = FW
      W = U
      FW = FU
      GO TO 20
   80 V = U
      FV = FU
      GO TO 20
C
C  end of main loop
C
   90 DFMIN = X
      RETURN
      END
      SUBROUTINE DFZERO(F,B,C,R,RE,AE,IFLAG) 

c*********************************************************************72
c
C***BEGIN PROLOGUE  DFZERO    
C***DATE WRITTEN   700901   (YYMMDD)    
C***REVISION DATE  890421   (YYMMDD)    
C***CATEGORY NO.  F1B    
C***KEYWORDS  BISECTION,DOUBLE PRECISION,NONLINEAR,ROOTS,ZEROS   
C***AUTHOR  SHAMPINE,L.F.,SNLA
C           WATTS,H.A.,SNLA   
C***PURPOSE  Search for a zero of a function F(X) in a given
C            interval (B,C).  It is designed primarily for problems   
C            where F(B) and F(C) have opposite signs.  
C***DESCRIPTION
C    
C       **** Double Precision version of FZERO ****    
C    
C     Based on a method by T J Dekker   
C     written by L F Shampine and H A Watts  
C    
C            DFZERO searches for a zero of a function F(X) between    
C            the given values B and C until the width of the interval 
C            (B,C) has collapsed to within a tolerance specified by   
C            the stopping criterion, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).    
C            The method used is an efficient combination of bisection 
C            and the secant rule.  
C    
C     Description Of Arguments
C    
C     F,B,C,R,RE and AE are DOUBLE PRECISION input parameters.   
C     B and C are DOUBLE PRECISION output parameters and IFLAG (flagged    
C        by an * below). 
C    
C        F     - Name of the DOUBLE PRECISION valued external function.    
C                This name must be in an EXTERNAL statement in the    
C                calling program.  F must be a function of one double 
C                precision argument.    
C    
C       *B     - One end of the interval (B,C).  The value returned for    
C                B usually is the better approximation to a zero of F.
C    
C       *C     - The other end of the interval (B,C)   
C    
C        R     - A (better) guess of a zero of F which could help in  
C                speeding up convergence.  If F(B) and F(R) have 
C                opposite signs, a root will be found in the interval 
C                (B,R); if not, but F(R) and F(C) have opposite  
C                signs, a root will be found in the interval (R,C);   
C                otherwise, the interval (B,C) will be searched for a 
C                possible root.  When no better guess is known, it is 
C                recommended that r be set to B or C; because if R is 
C                not interior to the interval (B,C), it will be ignored.   
C    
C        RE    - Relative error used for RW in the stopping criterion.
C                If the requested RE is less than machine precision,  
C                then RW is set to approximately machine precision.   
C    
C        AE    - Absolute error used in the stopping criterion.  If the    
C                given interval (B,C) contains the origin, then a
C                nonzero value should be chosen for AE.
C    
C       *IFLAG - A status code.  User must check IFLAG after each call.    
C                Control returns to the user from FZERO in all cases. 
C                XERROR does not process diagnostics in these cases.  
C    
C                1  B is within the requested tolerance of a zero.    
C                   The interval (B,C) collapsed to the requested
C                   tolerance, the function changes sign in (B,C), and
C                   F(X) decreased in magnitude as (B,C) collapsed.   
C    
C                2  F(B) = 0.  However, the interval (B,C) may not have    
C                   collapsed to the requested tolerance.   
C    
C                3  B may be near a singular point of F(X). 
C                   The interval (B,C) collapsed to the requested tol-
C                   erance and the function changes sign in (B,C), but
C                   F(X) increased in magnitude as (B,C) collapsed,i.e.    
C                     abs(F(B out)) .GT. max(abs(F(B in)),abs(F(C in)))    
C    
C                4  No change in sign of F(X) was found although the  
C                   interval (B,C) collapsed to the requested tolerance.   
C                   The user must examine this case and decide whether
C                   B is near a local minimum of F(X), or B is near a 
C                   zero of even multiplicity, or neither of these.   
C    
C                5  Too many (.GT. 500) function evaluations used.    
C***REFERENCES  L. F. SHAMPINE AND H. A. WATTS, *FZERO, A ROOT-SOLVING
C                 CODE*, SC-TM-70-631, SEPTEMBER 1970. 
C               T. J. DEKKER, *FINDING A ZERO BY MEANS OF SUCCESSIVE  
C                 LINEAR INTERPOLATION*, 'CONSTRUCTIVE ASPECTS OF THE 
C                 FUNDAMENTAL THEOREM OF ALGEBRA', EDITED BY B. DEJON 
C                 P. HENRICI, 1969.
C***ROUTINES CALLED  D1MACH   
C***END PROLOGUE  DFZERO 
C
      DOUBLE PRECISION A,ACBS,ACMB,AE,AW,B,C,CMB,D1MACH,ER,
     +                 F,FA,FB,FC,FX,FZ,P,Q,R,RE,RW,T,TOL,Z 
      INTEGER IC,IFLAG,KOUNT  
C    
C     ER IS TWO TIMES THE COMPUTER UNIT ROUNDOFF VALUE WHICH IS  
C     DEFINED HERE BY THE FUNCTION D1MACH.   
C    
C***FIRST EXECUTABLE STATEMENT  DFZERO   
      ER = 2.0D0 * D1MACH(4)  
C    
C     INITIALIZE    
C    
      Z = R 
      IF (R .LE. MIN(B,C)  .OR.  R .GE. MAX(B,C)) Z = C  
      RW = MAX(RE,ER)    
      AW = MAX(AE,0.D0)   
      IC = 0
      T = Z 
      FZ = F(T)  
      FC = FZ    
      T = B 
      FB = F(T)  
      KOUNT = 2  
      IF (SIGN(1.0D0,FZ) .EQ. SIGN(1.0D0,FB)) GO TO 1
      C = Z 
      GO TO 2  
    1 IF (Z .EQ. C) GO TO 2 
      T = C 
      FC = F(T)  
      KOUNT = 3  
      IF (SIGN(1.0D0,FZ) .EQ. SIGN(1.0D0,FC)) GO TO 2
      B = Z 
      FB = FZ    
    2 A = C 
      FA = FC    
      ACBS = ABS(B-C) 
      FX = MAX(ABS(FB),ABS(FC))    
C    
    3 IF (ABS(FC) .GE. ABS(FB)) GO TO 4 
C     PERFORM INTERCHANGE
      A = B 
      FA = FB    
      B = C 
      FB = FC    
      C = A 
      FC = FA    
C    
    4 CMB = 0.5D0*(C-B) 
      ACMB = ABS(CMB) 
      TOL = RW*ABS(B) + AW   
C    
C     TEST STOPPING CRITERION AND FUNCTION COUNT  
C    
      IF (ACMB .LE. TOL) GO TO 10  
      IF (FB .EQ. 0.D0) GO TO 11 
      IF (KOUNT .GE. 500) GO TO 14    
C    
C     CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q   
C     WHERE WE ARRANGE P .GE. 0.   
C     THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.   
C    
      P = (B-A)*FB    
      Q = FA - FB  
      IF (P .GE. 0.D0) GO TO 5  
      P = -P
      Q = -Q
C    
C     UPDATE A AND CHECK FOR SATISFACTORY REDUCTION    
C     IN THE SIZE OF THE BRACKETING INTERVAL.
C     IF NOT, PERFORM BISECTION.   
C    
    5 A = B 
      FA = FB    
      IC = IC + 1  
      IF (IC .LT. 4) GO TO 6  
      IF (8.0D0*ACMB .GE. ACBS) GO TO 8    
      IC = 0
      ACBS = ACMB
C    
C     TEST FOR TOO SMALL A CHANGE  
C    
    6 IF (P .GT. ABS(Q)*TOL) GO TO 7    
C    
C     INCREMENT BY TOLERANCE  
C    
      B = B + SIGN(TOL,CMB)  
      GO TO 9  
C    
C     ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.
C    
    7 IF (P .GE. CMB*Q) GO TO 8    
C    
C     USE SECANT RULE    
C    
      B = B + P/Q  
      GO TO 9  
C    
C     USE BISECTION   (C+B)/2 
C    
    8 B = B + CMB
C    
C     HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
C    
    9 T = B 
      FB = F(T)  
      KOUNT = KOUNT + 1 
C    
C     DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION 
C    
      IF (SIGN(1.0D0,FB) .NE. SIGN(1.0D0,FC)) GO TO 3 
      C = A 
      FC = FA    
      GO TO 3  
C    
C    
C     FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG 
C    
   10 IF (SIGN(1.0D0,FB) .EQ. SIGN(1.0D0,FC)) GO TO 13
      IF (ABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN   
   11 IFLAG = 2
      RETURN   
   12 IFLAG = 3
      RETURN   
   13 IFLAG = 4
      RETURN   
   14 IFLAG = 5
      RETURN   
      END 
      SUBROUTINE DGAMLM(XMIN,XMAX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGAMLM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C7A,R2
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(GAMLIM-S DGAMLM-D),
C             COMPLETE GAMMA FUNCTION,GAMMA FUNCTION,LIMITS,
C             SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the d.p. minimum and maximum bounds for X in
C            GAMMA(X).
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in gamma(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   double precision minimum legal value of X in gamma(X).  Any
C        smaller value of X might result in underflow.
C XMAX   double precision maximum legal value of X in gamma(X).  Any
C        larger value of X might cause overflow.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,XERROR
C***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
C***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = DLOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = DLOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (DABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMIN', 27, 1, 2)
C
 20   XMIN = -XMIN + 0.01D0
C
      ALNBIG = DLOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = DLOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (DABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERROR ( 'DGAMLM  UNABLE TO FIND XMAX', 27, 2, 2)
C
 40   XMAX = XMAX - 0.01D0
      XMIN = DMAX1 (XMIN, -XMAX+1.D0)

      RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMMA(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGAMMA
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  LIBRARY=SLATEC(FNLIB),
C             TYPE=DOUBLE PRECISION(GAMMA-S DGAMMA-D CGAMMA-C),
C             COMPLETE GAMMA FUNCTION,GAMMA FUNCTION,SPECIAL FUNCTIONS
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Compute the complete Gamma function.
C***DESCRIPTION
C
C DGAMMA(X) calculates the double precision complete Gamma function
C for double precision argument X.
C
C Series for GAM        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.00
C                                    decimal places required  32.05
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DCSEVL,DGAMLM,INITDS,XERROR
C***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX,
     1  XMIN, Y,  DINT, D9LGMC, DCSEVL, D1MACH
C
      SAVE GAM CS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL
      DATA GAMCS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAMCS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAMCS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAMCS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAMCS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAMCS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAMCS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAMCS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAMCS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAMCS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAMCS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAMCS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAMCS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAMCS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAMCS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAMCS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAMCS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAMCS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAMCS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAMCS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAMCS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAMCS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAMCS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAMCS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAMCS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAMCS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAMCS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAMCS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAMCS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAMCS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAMCS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAMCS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAMCS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAMCS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAMCS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAMCS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAMCS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAMCS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAMCS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAMCS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAMCS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAMCS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA NGAM, XMIN, XMAX, DXREL / 0, 0.D0,0.D0,0.D0 /
C***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (NGAM.NE.0) GO TO 10
      NGAM = INITDS (GAMCS, 42, 0.1*SNGL(D1MACH(3)) )
C
      CALL DGAMLM (XMIN, XMAX)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS(X)
      IF (Y.GT.10.D0) GO TO 50
C
C COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
C GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
C
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - DBLE(FLOAT(N))
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.0
C
      N = -N
      IF (X.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS 0', 14, 4, 2)
      IF (X.LT.0.0 .AND. X+DBLE(FLOAT(N-2)).EQ.0.D0) CALL XERROR ( 'DGAM
     1MA  X IS A NEGATIVE INTEGER', 31, 4, 2)
      IF (X.LT.(-0.5D0) .AND. DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL
     1  XERROR ( 'DGAMMA  ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NE
     2GATIVE INTEGER', 68, 1, 1)
C
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+DBLE(FLOAT(I-1)) )
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DGAMMA = (Y+DBLE(FLOAT(I))) * DGAMMA
 40   CONTINUE
      RETURN
C
C GAMMA(X) FOR DABS(X) .GT. 10.0.  RECALL Y = DABS(X).
C
 50   IF (X.GT.XMAX) CALL XERROR ( 'DGAMMA  X SO BIG GAMMA OVERFLOWS',
     1  32, 3, 2)
C
      DGAMMA = 0.D0
      IF (X.LT.XMIN) CALL XERROR ( 'DGAMMA  X SO SMALL GAMMA UNDERFLOWS'
     1  , 35, 2, 1)
      IF (X.LT.XMIN) RETURN
C
      DGAMMA = DEXP ((Y-0.5D0)*DLOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
C
      IF (DABS((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DGAMMA  ANS
     1WER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER'  , 61, 1, 1)
C
      SINPIY = DSIN (PI*Y)
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DGAMMA  X IS A NEGATIVE INTEGER
     1', 31, 4, 2)
C
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
C
      RETURN
      END
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGBFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  BANDED,DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,
C             MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision BAND matrix by elimination.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C     Fortran MAX0,MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      DOUBLE PRECISION ABD(LDA,1)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN0(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN0(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGBSL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A2
C***KEYWORDS  BANDED,DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,
C             SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision BAND system  A*X=B or
C            TRANS(A)*X=B using the factors computed by DGBCO or DGBFA.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C     Fortran MIN0
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      DOUBLE PRECISION ABD(LDA,1),B(1)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN0(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN0(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN0(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DGDRVD(NR,N,X,F,G,A,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,SC,WRK1,WRK2,WRK3,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE (XPLS) BY THE DOUBLE DOGLEG METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT  AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN
C                  IN LOWER TRIANGULAR PART AND DIAGONAL
C P(N)         --> NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C IRETCD      <--  RETURN CODE
C                    =0 SATISFACTORY XPLS FOUND
C                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
C                       DISTINCT FROM X
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C SC(N)        --> WORKSPACE [CURRENT STEP]
C WRK1(N)      --> WORKSPACE (AND PLACE HOLDING ARGUMENT TO TRGUPD)
C WRK2(N)      --> WORKSPACE
C WRK3(N)      --> WORKSPACE
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),XPLS(N),G(N),P(N)
      DIMENSION SX(N)
      DIMENSION SC(N),WRK1(N),WRK2(N),WRK3(N)
      DIMENSION A(NR,1)
      LOGICAL FSTDOG,NWTAKE,MXTAKE
      EXTERNAL FCN
C
      IRETCD=4
      FSTDOG=.TRUE.
      TMP=0.D0
      DO 5 I=1,N
        TMP=TMP+SX(I)*SX(I)*P(I)*P(I)
    5 CONTINUE
      RNWTLN=SQRT(TMP)
C$    WRITE(IPR,954) RNWTLN
C
  100 CONTINUE
C
C FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM
      CALL DGSTPD(NR,N,G,A,P,SX,RNWTLN,DLT,NWTAKE,FSTDOG,
     +     WRK1,WRK2,CLN,ETA,SC,IPR,STEPMX)
C
C CHECK NEW POINT AND UPDATE TRUST REGION
      CALL TRGUPD(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,DLT,
     +     IRETCD,WRK3,FPLSP,XPLS,FPLS,MXTAKE,IPR,2,WRK1)
      IF(IRETCD.LE.1) RETURN
      GO TO 100
  950 FORMAT(42H DGDRVD    INITIAL TRUST REGION NOT GIVEN.,
     +       22H  COMPUTE CAUCHY STEP.)
  951 FORMAT(18H DGDRVD    ALPHA =,E20.13/
     +       18H DGDRVD    BETA  =,E20.13/
     +       18H DGDRVD    DLT   =,E20.13/
     +       18H DGDRVD    NWTAKE=,L1    )
  952 FORMAT(28H DGDRVD    CURRENT STEP (SC))
  954 FORMAT(18H0DGDRVD    RNWTLN=,E20.13)
  955 FORMAT(14H DGDRVD       ,5(E20.13,3X))
      END
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)                           

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGECO                                              
C***DATE WRITTEN   780814   (YYMMDD)                                   
C***REVISION DATE  861211   (YYMMDD)                                   
C***CATEGORY NO.  D2A1                                                 
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),                                 
C             TYPE=DOUBLE PRECISION(SGECO-S DGECO-D CGECO-C),          
C             CONDITION NUMBER,GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,   
C             MATRIX FACTORIZATION                                     
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                           
C***PURPOSE  Factors a double precision matrix by Gaussian elimination 
C            and estimates the condition of the matrix.                
C***DESCRIPTION                                                        
C                                                                      
C     DGECO factors a double precision matrix by Gaussian elimination  
C     and estimates the condition of the matrix.                       
C                                                                      
C     If  RCOND  is not needed, DGEFA is slightly faster.              
C     To solve  A*X = B , follow DGECO by DGESL.                       
C     To compute  INVERSE(A)*C , follow DGECO by DGESL.                
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.              
C     To compute  INVERSE(A) , follow DGECO by DGEDI.                  
C                                                                      
C     On Entry                                                         
C                                                                      
C        A       DOUBLE PRECISION(LDA, N)                              
C                the matrix to be factored.                            
C                                                                      
C        LDA     INTEGER                                               
C                the leading dimension of the array  A .               
C                                                                      
C        N       INTEGER                                               
C                the order of the matrix  A .                          
C                                                                      
C     On Return                                                        
C                                                                      
C        A       an upper triangular matrix and the multipliers        
C                which were used to obtain it.                         
C                The factorization can be written  A = L*U  where      
C                L  is a product of permutation and unit lower         
C                triangular matrices and  U  is upper triangular.      
C                                                                      
C        IPVT    INTEGER(N)                                            
C                an INTEGER vector of pivot indices.                   
C                                                                      
C        RCOND   DOUBLE PRECISION                                      
C                an estimate of the reciprocal condition of  A .       
C                For the system  A*X = B , relative perturbations      
C                in  A  and  B  of size  EPSILON  may cause            
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression    
C                           1.0 + RCOND .EQ. 1.0                       
C                is true, then  A  may be singular to working          
C                precision.  In particular,  RCOND  is zero  if        
C                exact singularity is detected or the estimate         
C                underflows.                                           
C                                                                      
C        Z       DOUBLE PRECISION(N)                                   
C                a work vector whose contents are usually unimportant. 
C                If  A  is close to a singular matrix, then  Z  is     
C                an approximate null vector in the sense that          
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                   
C                                                                      
C     LINPACK.  This version dated 08/14/78 .                          
C     Cleve Moler, University of New Mexico, Argonne National Lab.     
C                                                                      
C     Subroutines and Functions                                        
C                                                                      
C     LINPACK DGEFA                                                    
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                      
C     Fortran DABS,DMAX1,DSIGN                                         
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,   
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                  
C***ROUTINES CALLED  DASUM,DAXPY,DDOT,DGEFA,DSCAL                      
C***END PROLOGUE  DGECO                                                
      INTEGER LDA,N,IPVT(1)                                            
      DOUBLE PRECISION A(LDA,1),Z(1)                                   
      DOUBLE PRECISION RCOND                                           
C                                                                      
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                          
      INTEGER INFO,J,K,KB,KP1,L                                        
C                                                                      
C     COMPUTE 1-NORM OF A                                              
C                                                                      
C***FIRST EXECUTABLE STATEMENT  DGECO                                  
      ANORM = 0.0D0                                                    
      DO 10 J = 1, N                                                   
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))                        
   10 CONTINUE                                                         
C                                                                      
C     FACTOR                                                           
C                                                                      
      CALL DGEFA(A,LDA,N,IPVT,INFO)                                    
C                                                                      
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .             
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E . 
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE     
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE 
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID   
C     OVERFLOW.                                                        
C                                                                      
C     SOLVE TRANS(U)*W = E                                             
C                                                                      
      EK = 1.0D0                                                       
      DO 20 J = 1, N                                                   
         Z(J) = 0.0D0                                                  
   20 CONTINUE                                                         
      DO 100 K = 1, N                                                  
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                     
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30                 
            S = DABS(A(K,K))/DABS(EK-Z(K))                             
            CALL DSCAL(N,S,Z,1)                                        
            EK = S*EK                                                  
   30    CONTINUE                                                      
         WK = EK - Z(K)                                                
         WKM = -EK - Z(K)                                              
         S = DABS(WK)                                                  
         SM = DABS(WKM)                                                
         IF (A(K,K) .EQ. 0.0D0) GO TO 40                               
            WK = WK/A(K,K)                                             
            WKM = WKM/A(K,K)                                           
         GO TO 50                                                      
   40    CONTINUE                                                      
            WK = 1.0D0                                                 
            WKM = 1.0D0                                                
   50    CONTINUE                                                      
         KP1 = K + 1                                                   
         IF (KP1 .GT. N) GO TO 90                                      
            DO 60 J = KP1, N                                           
               SM = SM + DABS(Z(J)+WKM*A(K,J))                         
               Z(J) = Z(J) + WK*A(K,J)                                 
               S = S + DABS(Z(J))                                      
   60       CONTINUE                                                   
            IF (S .GE. SM) GO TO 80                                    
               T = WKM - WK                                            
               WK = WKM                                                
               DO 70 J = KP1, N                                        
                  Z(J) = Z(J) + T*A(K,J)                               
   70          CONTINUE                                                
   80       CONTINUE                                                   
   90    CONTINUE                                                      
         Z(K) = WK                                                     
  100 CONTINUE                                                         
      S = 1.0D0/DASUM(N,Z,1)                                           
      CALL DSCAL(N,S,Z,1)                                              
C                                                                      
C     SOLVE TRANS(L)*Y = W                                             
C                                                                      
      DO 120 KB = 1, N                                                 
         K = N + 1 - KB                                                
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)     
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                          
            S = 1.0D0/DABS(Z(K))                                       
            CALL DSCAL(N,S,Z,1)                                        
  110    CONTINUE                                                      
         L = IPVT(K)                                                   
         T = Z(L)                                                      
         Z(L) = Z(K)                                                   
         Z(K) = T                                                      
  120 CONTINUE                                                         
      S = 1.0D0/DASUM(N,Z,1)                                           
      CALL DSCAL(N,S,Z,1)                                              
C                                                                      
      YNORM = 1.0D0                                                    
C                                                                      
C     SOLVE L*V = Y                                                    
C                                                                      
      DO 140 K = 1, N                                                  
         L = IPVT(K)                                                   
         T = Z(L)                                                      
         Z(L) = Z(K)                                                   
         Z(K) = T                                                      
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)           
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                          
            S = 1.0D0/DABS(Z(K))                                       
            CALL DSCAL(N,S,Z,1)                                        
            YNORM = S*YNORM                                            
  130    CONTINUE                                                      
  140 CONTINUE                                                         
      S = 1.0D0/DASUM(N,Z,1)                                           
      CALL DSCAL(N,S,Z,1)                                              
      YNORM = S*YNORM                                                  
C                                                                      
C     SOLVE  U*Z = V                                                   
C                                                                      
      DO 160 KB = 1, N                                                 
         K = N + 1 - KB                                                
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150                   
            S = DABS(A(K,K))/DABS(Z(K))                                
            CALL DSCAL(N,S,Z,1)                                        
            YNORM = S*YNORM                                            
  150    CONTINUE                                                      
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                     
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                           
         T = -Z(K)                                                     
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                             
  160 CONTINUE                                                         
C     MAKE ZNORM = 1.0                                                 
      S = 1.0D0/DASUM(N,Z,1)                                           
      CALL DSCAL(N,S,Z,1)                                              
      YNORM = S*YNORM                                                  
C                                                                      
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                        
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                              
      RETURN                                                           
      END                                                              
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)   

c*********************************************************************72
c                               
C***BEGIN PROLOGUE  DGEFA                                                  
C***DATE WRITTEN   780814   (YYMMDD)                                       
C***REVISION DATE  861211   (YYMMDD)                                       
C***CATEGORY NO.  D2A1                                                     
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),                                     
C             TYPE=DOUBLE PRECISION(SGEFA-S DGEFA-D CGEFA-C),              
C             GENERAL MATRIX,LINEAR ALGEBRA,MATRIX,MATRIX FACTORIZATION    
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                               
C***PURPOSE  Factors a double precision matrix by Gaussian elimination.    
C***DESCRIPTION                                                            
C                                                                          
C     DGEFA factors a double precision matrix by Gaussian elimination.     
C                                                                          
C     DGEFA is usually called by DGECO, but it can be called               
C     directly with a saving in time if  RCOND  is not needed.             
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .                      
C                                                                          
C     On Entry                                                             
C                                                                          
C        A       DOUBLE PRECISION(LDA, N)                                  
C                the matrix to be factored.                                
C                                                                          
C        LDA     INTEGER                                                   
C                the leading dimension of the array  A .                   
C                                                                          
C        N       INTEGER                                                   
C                the order of the matrix  A .                              
C                                                                          
C     On Return                                                            
C                                                                          
C        A       an upper triangular matrix and the multipliers            
C                which were used to obtain it.                             
C                The factorization can be written  A = L*U  where          
C                L  is a product of permutation and unit lower             
C                triangular matrices and  U  is upper triangular.          
C                                                                          
C        IPVT    INTEGER(N)                                                
C                an integer vector of pivot indices.                       
C                                                                          
C        INFO    INTEGER                                                   
C                = 0  normal value.                                        
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error          
C                     condition for this subroutine, but it does           
C                     indicate that DGESL or DGEDI will divide by zero     
C                     if called.  Use  RCOND  in DGECO for a reliable      
C                     indication of singularity.                           
C                                                                          
C     LINPACK.  This version dated 08/14/78 .                              
C     Cleve Moler, University of New Mexico, Argonne National Lab.         
C                                                                          
C     Subroutines and Functions                                            
C                                                                          
C     BLAS DAXPY,DSCAL,IDAMAX                                              
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,       
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                      
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX                                    
C***END PROLOGUE  DGEFA                                                    
      INTEGER LDA,N,IPVT(1),INFO                                           
      DOUBLE PRECISION A(LDA,1)                                            
C                                                                          
      DOUBLE PRECISION T                                                   
      INTEGER IDAMAX,J,K,KP1,L,NM1                                         
C                                                                          
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                           
C                                                                          
C***FIRST EXECUTABLE STATEMENT  DGEFA                                      
      INFO = 0                                                             
      NM1 = N - 1                                                          
      IF (NM1 .LT. 1) GO TO 70                                             
      DO 60 K = 1, NM1                                                     
         KP1 = K + 1                                                       
C                                                                          
C        FIND L = PIVOT INDEX                                              
C                                                                          
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                                
         IPVT(K) = L                                                       
C                                                                          
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED             
C                                                                          
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                   
C                                                                          
C           INTERCHANGE IF NECESSARY                                       
C                                                                          
            IF (L .EQ. K) GO TO 10                                         
               T = A(L,K)                                                  
               A(L,K) = A(K,K)                                             
               A(K,K) = T                                                  
   10       CONTINUE                                                       
C                                                                          
C           COMPUTE MULTIPLIERS                                            
C                                                                          
            T = -1.0D0/A(K,K)                                              
            CALL DSCAL(N-K,T,A(K+1,K),1)                                   
C                                                                          
C           ROW ELIMINATION WITH COLUMN INDEXING                           
C                                                                          
            DO 30 J = KP1, N                                               
               T = A(L,J)                                                  
               IF (L .EQ. K) GO TO 20                                      
                  A(L,J) = A(K,J)                                          
                  A(K,J) = T                                               
   20          CONTINUE                                                    
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                     
   30       CONTINUE                                                       
         GO TO 50                                                          
   40    CONTINUE                                                          
            INFO = K                                                       
   50    CONTINUE                                                          
   60 CONTINUE                                                             
   70 CONTINUE                                                             
      IPVT(N) = N                                                          
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                      
      RETURN                                                               
      END                                                                  
      SUBROUTINE DGEFS(A,LDA,N,V,ITASK,IND,WORK,IWORK,RCOND)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DGEFS
C***DATE WRITTEN   800326   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D2A1
C***AUTHOR  VOORHEES, E., (LANL)
C***PURPOSE  DGEFS solves a GENERAL double precision
C            NXN system of linear equations.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by D. Kahaner, C. Moler, S. Nash
C          Prentice Hall 1988
C
C    Subroutine DGEFS solves a general NxN system of double
C    precision linear equations using LINPACK subroutines DGECO
C    and DGESL.  That is, if A is an NxN double precision matrix
C    and if X and B are double precision N-vectors, then DGEFS
C    solves the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to find the
C    solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK.GT.1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N and IWORK must not have been altered by the user follow-
C    ing factorization (ITASK=1).  IND will not be changed by DGEFS
C    in this case.
C
C  Argument Description ***
C
C    A      DOUBLE PRECISION(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.
C             on return, an upper triangular matrix U and the
C               multipliers necessary to construct a matrix L
C               so that A=L*U.
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  The first N elements of
C             the array A are the elements of the first column of
C             the matrix A.  N must be greater than or equal to 1.
C             (terminal error message IND=-2)
C    V      DOUBLE PRECISION(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             If ITASK=1, the matrix A is factored and then the
C               linear equation is solved.
C             If ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A and IWORK.
C             If ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   DOUBLE PRECISION(N)
C             a singly subscripted array of dimension at least N.
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
C
C               Note-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERRWV in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERROR.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C***REFERENCES  SUBROUTINE DGEFS WAS DEVELOPED BY GROUP C-3, LOS ALAMOS
C                 SCIENTIFIC LABORATORY, LOS ALAMOS, NM 87545.
C                 THE LINPACK SUBROUTINES USED BY DGEFS ARE DESCRIBED IN
C                 DETAIL IN THE *LINPACK USERS GUIDE* PUBLISHED BY
C                 THE SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS
C                 (SIAM) DATED 1979.
C***ROUTINES CALLED  D1MACH,DGECO,DGESL,XERROR,XERRWV
C***END PROLOGUE  DGEFS
C
      INTEGER LDA,N,ITASK,IND,IWORK(N)
      DOUBLE PRECISION A(LDA,N),V(N),WORK(N),D1MACH
      DOUBLE PRECISION RCOND
C***FIRST EXECUTABLE STATEMENT  DGEFS
      IF (LDA.LT.N)  GO TO 101
      IF (N.LE.0)  GO TO 102
      IF (ITASK.LT.1) GO TO 103
      IF (ITASK.GT.1) GO TO 20
C
C     FACTOR MATRIX A INTO LU
      CALL DGECO(A,LDA,N,IWORK,RCOND,WORK)
C
C     CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
      IF (RCOND.EQ.0.0D0)  GO TO 104
C
C     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
      IND=-IDINT(DLOG10(D1MACH(4)/RCOND))
C
C     CHECK FOR IND GREATER THAN ZERO
      IF (IND.GT.0)  GO TO 20
      IND=-10
      CALL XERROR( 'DGEFS ERROR (IND=-10) -- SOLUTION MAY HAVE NO SIGNIF
     1ICANCE',58,-10,0)
C
C     SOLVE AFTER FACTORING
   20 CALL DGESL(A,LDA,N,IWORK,V,0)
      RETURN
C
C     IF LDA.LT.N, IND=-1, TERMINAL XERRWV MESSAGE
  101 IND=-1
      CALL XERRWV( 'DGEFS ERROR (IND=-1) -- LDA=I1 IS LESS THAN N=I2',
     148,-1,1,2,LDA,N,0,0,0)
      RETURN
C
C     IF N.LT.1, IND=-2, TERMINAL XERRWV MESSAGE
  102 IND=-2
      CALL XERRWV( 'DGEFS ERROR (IND=-2) -- N=I1 IS LESS THAN 1',
     143,-2,1,1,N,0,0,0,0)
      RETURN
C
C     IF ITASK.LT.1, IND=-3, TERMINAL XERRWV MESSAGE
  103 IND=-3
      CALL XERRWV( 'DGEFS ERROR (IND=-3) -- ITASK=I1 IS LESS THAN 1',
     147,-3,1,1,ITASK,0,0,0,0)
      RETURN
C
C     IF SINGULAR MATRIX, IND=-4, TERMINAL XERRWV MESSAGE
  104 IND=-4
      CALL XERRWV( 'DGEFS ERROR (IND=-4) -- SINGULAR MATRIX A - NO SOLUT
     1ION',55,-4,1,0,0,0,0,0,0)
      RETURN
C
      END
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)

c*********************************************************************72
c                           
C***BEGIN PROLOGUE  DGESL                                               
C***DATE WRITTEN   780814   (YYMMDD)                                    
C***REVISION DATE  861211   (YYMMDD)                                    
C***CATEGORY NO.  D2A1                                                  
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),                                  
C             TYPE=DOUBLE PRECISION(SGESL-S DGESL-D CGESL-C),           
C             LINEAR ALGEBRA,MATRIX,SOLVE                               
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)                            
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B 
C            using the factors computed by DGECO or DGEFA.              
C***DESCRIPTION                                                         
C                                                                       
C     DGESL solves the double precision system                          
C     A * X = B  or  TRANS(A) * X = B                                   
C     using the factors computed by DGECO or DGEFA.                     
C                                                                       
C     On Entry                                                          
C                                                                       
C        A       DOUBLE PRECISION(LDA, N)                               
C                the output from DGECO or DGEFA.                        
C                                                                       
C        LDA     INTEGER                                                

C                the leading dimension of the array  A .                
C                                                                       
C        N       INTEGER                                                
C                the order of the matrix  A .                           
C                                                                       
C        IPVT    INTEGER(N)                                             
C                the pivot vector from DGECO or DGEFA.                  
C                                                                       
C        B       DOUBLE PRECISION(N)                                    
C                the right hand side vector.                            
C                                                                       
C        JOB     INTEGER                                                
C                = 0         to solve  A*X = B ,                        
C                = nonzero   to solve  TRANS(A)*X = B  where            
C                            TRANS(A)  is the transpose.                
C                                                                       
C     On Return                                                         
C                                                                       
C        B       the solution vector  X .                               
C                                                                       
C     Error Condition                                                   
C                                                                       
C        A division by zero will occur if the input factor contains a   
C        zero on the diagonal.  Technically this indicates singularity  
C        but it is often caused by improper arguments or improper       
C        setting of LDA .  It will not occur if the subroutines are     
C        called correctly and if DGECO has set RCOND .GT. 0.0           
C        or DGEFA has set INFO .EQ. 0 .                                 
C                                                                       
C     To compute  INVERSE(A) * C  where  C  is a matrix                 
C     with  P  columns                                                  
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            
C           IF (RCOND is too small) GO TO ...                           
C           DO 10 J = 1, P                                              
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        
C        10 CONTINUE                                                    
C                                                                       
C     LINPACK.  This version dated 08/14/78 .                           
C     Cleve Moler, University of New Mexico, Argonne National Lab.      
C                                                                       
C     Subroutines and Functions                                         
C                                                                       
C     BLAS DAXPY,DDOT                                                   
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   
C***ROUTINES CALLED  DAXPY,DDOT                                         
C***END PROLOGUE  DGESL                                                 
      INTEGER LDA,N,IPVT(1),JOB                                         
      DOUBLE PRECISION A(LDA,1),B(1)                                    
C                                                                       
      DOUBLE PRECISION DDOT,T                                           
      INTEGER K,KB,L,NM1                                                
C***FIRST EXECUTABLE STATEMENT  DGESL                                   
      NM1 = N - 1                                                       
      IF (JOB .NE. 0) GO TO 50                                          
C                                                                       
C        JOB = 0 , SOLVE  A * X = B                                     
C        FIRST SOLVE  L*Y = B                                           
C                                                                       
         IF (NM1 .LT. 1) GO TO 30                                       
         DO 20 K = 1, NM1                                               
            L = IPVT(K)                                                 
            T = B(L)                                                    
            IF (L .EQ. K) GO TO 10                                      
               B(L) = B(K)                                              
               B(K) = T                                                 
   10       CONTINUE                                                    
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)                       
   20    CONTINUE                                                       
   30    CONTINUE                                                       
C                                                                       
C        NOW SOLVE  U*X = Y                                             
C                                                                       
         DO 40 KB = 1, N                                                
            K = N + 1 - KB                                              
            B(K) = B(K)/A(K,K)                                          
            T = -B(K)                                                   
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                           
   40    CONTINUE                                                       
      GO TO 100                                                         
   50 CONTINUE                                                          
C                                                                       
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         
C        FIRST SOLVE  TRANS(U)*Y = B                                    
C                                                                       
         DO 60 K = 1, N                                                 
            T = DDOT(K-1,A(1,K),1,B(1),1)                               
            B(K) = (B(K) - T)/A(K,K)                                    
   60    CONTINUE                                                       
C                                                                       
C        NOW SOLVE TRANS(L)*X = Y                                       
C                                                                       
         IF (NM1 .LT. 1) GO TO 90                                       
         DO 80 KB = 1, NM1                                              
            K = N - KB                                                  
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)                 
            L = IPVT(K)                                                 
            IF (L .EQ. K) GO TO 70                                      
               T = B(L)                                                 
               B(L) = B(K)                                              
               B(K) = T                                                 
   70       CONTINUE                                                    
   80    CONTINUE                                                       
   90    CONTINUE                                                       
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DGL15T(F,A,B,XL,XR,R,AE,RA,
     1                 RASC,FMIN,FMAX)

c*********************************************************************72
c
C
C***AUTHORS          ROBERT PIESSENS AND ELISE DE DONCKER
C                    APPL. MATH. AND PROGR. DIV.- K.U.LEUVEN
C                    DAVID KAHANER, NBS WASHINGTON
C
C***PURPOSE
C              TO COMPUTE I = INTEGRAL OF G(X) OVER (A,B),
C                             WITH ERROR ESTIMATE
C                         J = INTEGRAL OF ABS(G) OVER (A,B)
C              DOUBLE PRECISION VERSION OF GL15T
C
C***DESCRIPTION
C           PARAMETERS
C            ON ENTRY
C              F      - DOUBLE PRECISION
C                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
C                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS
C                       TO BE DECLARED E X T E R N A L IN THE
C                       CALLING PROGRAM.
C                       THE FUNCTION G(X) IS DEFINED TO BE
C                       G(X)=F(PHI(X))*PHIP(X)
C                       WHERE PHI(X) IS THE CUBIC GIVEN BY
C                       THE ARITHMETIC STATEMENT FUNCTION BELOW.
C                       PHIP(X) IS ITS DERIVATIVE.  THE VARIABLES
C                       XL AND XR ARE THE LEFT AND RIGHT ENDPOINTS
C                       OF A PARENT INTERVAL OF WHICH (A,B) IS A PART.
C
C              A      - DOUBLE PRECISION
C                       LOWER LIMIT OF INTEGRATION
C
C              B      - DOUBLE PRECISION
C                       UPPER LIMIT OF INTEGRATION
C
C              XL     - DOUBLE PRECISION
C              XR     - DOUBLE PRECISION
C                       LOWER AND UPPER LIMITS OF PARENT INTERVAL
C                       OF WHICH [A,B] IS A PART.
C
C            ON RETURN
C              R - DOUBLE PRECISION
C                       APPROXIMATION TO THE INTEGRAL I
C                       R IS COMPUTED BY APPLYING THE 15-POINT
C                       KRONROD RULE (RESK) OBTAINED BY OPTIMAL
C                       ADDITION OF ABSCISSAE TO THE 7-POINT GAUSS
C                       RULE (RESG).
C
C              AE - DOUBLE PRECISION
C                       ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
C                       WHICH SHOULD NOT EXCEED ABS(I-R)
C
C              RA - DOUBLE PRECISION
C                       APPROXIMATION TO THE INTEGRAL J
C
C              RASC - DOUBLE PRECISION
C                       APPROXIMATION TO THE INTEGRAL OF ABS(G-I/(B-A))
C                       OVER (A,B)
C
C              FMAX, FMIN - DOUBLE PRECISION
C                       MAX AND MIN VALUES OF THE FUNCTION F ON (A,B)
C           SUBROUTINES OR FUNCTIONS NEEDED
C                 - F (USER-PROVIDED FUNCTION)
C                 - DOUBLE PRECISION
C                 - FORTRAN ABS, MAX, MIN
C
C     .................................................................
C***END PROLOGUE  DGL15T
C
      SAVE EPMACH,UFLOW
      DOUBLE PRECISION A,AE,B,DHLGTH,EPMACH,F,FC,FMAX,FMIN,FSUM,FV1,FV2,
     *  FVAL1,FVAL2,
     *  HLGTH,PHI,PHIP,PHIU,R,D1MACH,RA,RASC,RESG,RESK,RESKH,SL,SR,
     *  UFLOW,WG,WGK,XGK
      DOUBLE PRECISION XL,XR,CENTR,ABSC,U
      INTEGER J,JTW,JTWM1
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1)
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
      DATA WG  (  1) / 0.1294849661 6886969327 0611432679 082 D0 /
      DATA WG  (  2) / 0.2797053914 8927666790 1467771423 780 D0 /
      DATA WG  (  3) / 0.3818300505 0511894495 0369775488 975 D0 /
      DATA WG  (  4) / 0.4179591836 7346938775 5102040816 327 D0 /
C
      DATA XGK (  1) / 0.9914553711 2081263920 6854697526 329 D0 /
      DATA XGK (  2) / 0.9491079123 4275852452 6189684047 851 D0 /
      DATA XGK (  3) / 0.8648644233 5976907278 9712788640 926 D0 /
      DATA XGK (  4) / 0.7415311855 9939443986 3864773280 788 D0 /
      DATA XGK (  5) / 0.5860872354 6769113029 4144838258 730 D0 /
      DATA XGK (  6) / 0.4058451513 7739716690 6606412076 961 D0 /
      DATA XGK (  7) / 0.2077849550 0789846760 0689403773 245 D0 /
      DATA XGK (  8) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0229353220 1052922496 3732008058 970 D0 /
      DATA WGK (  2) / 0.0630920926 2997855329 0700663189 204 D0 /
      DATA WGK (  3) / 0.1047900103 2225018383 9876322541 518 D0 /
      DATA WGK (  4) / 0.1406532597 1552591874 5189590510 238 D0 /
      DATA WGK (  5) / 0.1690047266 3926790282 6583426598 550 D0 /
      DATA WGK (  6) / 0.1903505780 6478540991 3256402421 014 D0 /
      DATA WGK (  7) / 0.2044329400 7529889241 4161999234 649 D0 /
      DATA WGK (  8) / 0.2094821410 8472782801 2999174891 714 D0 /
C
C
      PHI(U)=XR-(XR-XL)*U*U*(2.*U+3.)
      PHIP(U)=-6.*U*(U+1.)
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - R OF THE 7-POINT GAUSS FORMULA
C           RESK   - R OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT
      DATA EPMACH,UFLOW/0.0,0.0/
      IF(EPMACH.EQ.0.0) THEN
         EPMACH=D1MACH(4)
         UFLOW=D1MACH(1)
      ENDIF
C
      IF(XL.LT.XR)THEN
         SL=(XL)
         SR=(XR)
        ELSE
         SL=(XR)
         SR=(XL)
      ENDIF
      HLGTH = 0.5D+00*(B-A)
      CENTR = A+HLGTH
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      U=(CENTR-XR)/(XR-XL)
      PHIU=PHI(U)
      IF(PHIU.LE.SL .OR. PHIU.GE.SR) PHIU=CENTR
      FMIN=F(PHIU)
      FMAX=FMIN
      FC=FMIN*PHIP(U)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RA = ABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        U=(CENTR-ABSC-XR)/(XR-XL)
        PHIU=PHI(U)
        IF(PHIU.LE.SL .OR. PHIU.GE.SR) PHIU=CENTR
        FVAL1=F(PHIU)
        FMAX=MAX(FMAX,FVAL1)
        FMIN=MIN(FMIN,FVAL1)
        FVAL1=FVAL1*PHIP(U)
        U=(CENTR+ABSC-XR)/(XR-XL)
        PHIU=PHI(U)
        IF(PHIU.LE.SL .OR. PHIU.GE.SR) PHIU=CENTR
        FVAL2=F(PHIU)
        FMAX=MAX(FMAX,FVAL2)
        FMIN=MIN(FMIN,FVAL2)
        FVAL2=FVAL2*PHIP(U)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RA = RA+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        U=(CENTR-ABSC-XR)/(XR-XL)
        PHIU=PHI(U)
        IF(PHIU.LE.SL .OR. PHIU.GE.SR) PHIU=CENTR
        FVAL1=F(PHIU)
        FMAX=MAX(FMAX,FVAL1)
        FMIN=MIN(FMIN,FVAL1)
        FVAL1=FVAL1*PHIP(U)
        U=(CENTR+ABSC-XR)/(XR-XL)
        PHIU=PHI(U)
        IF(PHIU.LE.SL .OR. PHIU.GE.SR) PHIU=CENTR
        FVAL2=F(PHIU)
        FMAX=MAX(FMAX,FVAL2)
        FMIN=MIN(FMIN,FVAL2)
        FVAL2=FVAL2*PHIP(U)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RA = RA+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RASC = RASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      R = RESK*HLGTH
      RA = RA*DHLGTH
      RASC = RASC*DHLGTH
      AE = ABS((RESK-RESG)*HLGTH)
      IF(RASC.NE.0.0D+00.AND.AE.NE.0.0D+00)
     *  AE = RASC*MIN(0.1D+01,
     *  (0.2D+03*AE/RASC)**1.5D+00)
      IF(RA.GT.UFLOW/(0.5D+02*EPMACH)) AE = MAX(
     *  (EPMACH*0.5D+02)*RA,AE)
      RETURN
      END
      SUBROUTINE DGSTPD(NR,N,G,A,P,SX,RNWTLN,DLT,NWTAKE,FSTDOG,
     +     SSD,V,CLN,ETA,SC,IPR,STEPMX)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER PART AND DIAGONAL
C P(N)         --> NEWTON STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNWTLN       --> NEWTON STEP LENGTH
C DLT         <--> TRUST REGION RADIUS
C NWTAKE      <--> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C FSTDOG      <--> BOOLEAN, =.TRUE. IF ON FIRST LEG OF DOGLEG
C SSD(N)      <--> WORKSPACE [CAUCHY STEP TO THE MINIMUM OF THE
C                  QUADRATIC MODEL IN THE SCALED STEEPEST DESCENT
C                  DIRECTION] [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C V(N)        <--> WORKSPACE  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C CLN         <--> CAUCHY LENGTH
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C ETA              [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C SC(N)       <--  CURRENT STEP
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C
C INTERNAL VARIABLES
C ------------------
C CLN              LENGTH OF CAUCHY STEP
C
      DIMENSION G(N),P(N)
      DIMENSION SX(N)
      DIMENSION SC(N),SSD(N),V(N)
      DIMENSION A(NR,1)
      LOGICAL NWTAKE,FSTDOG
      IPR=IPR
C
C CAN WE TAKE NEWTON STEP
C
      IF(RNWTLN.GT.DLT) GO TO 100
C     IF(RNWTLN.LE.DLT)
C     THEN
        NWTAKE=.TRUE.
        DO 10 I=1,N
          SC(I)=P(I)
   10   CONTINUE
        DLT=RNWTLN
C$      WRITE(IPR,951)
        GO TO 700
C     ELSE
C
C NEWTON STEP TOO LONG
C CAUCHY STEP IS ON DOUBLE DOGLEG CURVE
C
  100   NWTAKE=.FALSE.
        IF(.NOT.FSTDOG) GO TO 200
C       IF(FSTDOG)
C       THEN
C
C         CALCULATE DOUBLE DOGLEG CURVE (SSD)
          FSTDOG=.FALSE.
          ALPHA=0.D0
          DO 110 I=1,N
            ALPHA=ALPHA + (G(I)*G(I))/(SX(I)*SX(I))
  110     CONTINUE
          BETA=0.D0
          DO 130 I=1,N
            TMP=0.D0
            DO 120 J=I,N
              TMP=TMP + (A(J,I)*G(J))/(SX(J)*SX(J))
  120       CONTINUE
            BETA=BETA+TMP*TMP
  130     CONTINUE
          DO 140 I=1,N
            SSD(I)=-(ALPHA/BETA)*G(I)/SX(I)
  140     CONTINUE
          CLN=ALPHA*SQRT(ALPHA)/BETA
          ETA=.2D0 + (.8D0*ALPHA*ALPHA)/(-BETA*DDOT(N,G,1,P,1))
          DO 150 I=1,N
            V(I)=ETA*SX(I)*P(I) - SSD(I)
  150     CONTINUE
          IF (DLT .EQ. (-1.0D0)) DLT = MIN(CLN, STEPMX)
C$        WRITE(IPR,954) ALPHA,BETA,CLN,ETA
C$        WRITE(IPR,955)
C$        WRITE(IPR,960) (SSD(I),I=1,N)
C$        WRITE(IPR,956)
C$        WRITE(IPR,960) (V(I),I=1,N)
C       ENDIF
  200   IF(ETA*RNWTLN.GT.DLT) GO TO 220
C       IF(ETA*RNWTLN .LE. DLT)
C       THEN
C
C         TAKE PARTIAL STEP IN NEWTON DIRECTION
C
          DO 210 I=1,N
            SC(I)=(DLT/RNWTLN)*P(I)
  210     CONTINUE
C$        WRITE(IPR,957)
          GO TO 700
C       ELSE
  220     IF(CLN.LT.DLT) GO TO 240
C         IF(CLN.GE.DLT)
C         THEN
C           TAKE STEP IN STEEPEST DESCENT DIRECTION
C
            DO 230 I=1,N
              SC(I)=(DLT/CLN)*SSD(I)/SX(I)
  230       CONTINUE
C$          WRITE(IPR,958)
            GO TO 700
C         ELSE
C
C           CALCULATE CONVEX COMBINATION OF SSD AND ETA*P
C           WHICH HAS SCALED LENGTH DLT
C
  240       DOT1=DDOT(N,V,1,SSD,1)
            DOT2=DDOT(N,V,1,V,1)
            ALAM=(-DOT1+SQRT((DOT1*DOT1)-DOT2*(CLN*CLN-DLT*DLT)))/DOT2
            DO 250 I=1,N
              SC(I)=(SSD(I) + ALAM*V(I))/SX(I)
  250       CONTINUE
C$          WRITE(IPR,959)
C         ENDIF
C       ENDIF
C     ENDIF
  700 CONTINUE
C$    WRITE(IPR,952) FSTDOG,NWTAKE,RNWTLN,DLT
C$    WRITE(IPR,953)
C$    WRITE(IPR,960) (SC(I),I=1,N)
      RETURN
C
  951 FORMAT(27H0DGSTPD    TAKE NEWTON STEP)
  952 FORMAT(18H DGSTPD    FSTDOG=,L1/
     +       18H DGSTPD    NWTAKE=,L1/
     +       18H DGSTPD    RNWTLN=,E20.13/
     +       18H DGSTPD    DLT   =,E20.13)
  953 FORMAT(28H DGSTPD    CURRENT STEP (SC))
  954 FORMAT(18H DGSTPD    ALPHA =,E20.13/
     +       18H DGSTPD    BETA  =,E20.13/
     +       18H DGSTPD    CLN   =,E20.13/
     +       18H DGSTPD    ETA   =,E20.13)
  955 FORMAT(28H DGSTPD    CAUCHY STEP (SSD))
  956 FORMAT(12H DGSTPD    V)
  957 FORMAT(48H0DGSTPD    TAKE PARTIAL STEP IN NEWTON DIRECTION)
  958 FORMAT(50H0DGSTPD    TAKE STEP IN STEEPEST DESCENT DIRECTION)
  959 FORMAT(39H0DGSTPD    TAKE CONVEX COMBINATION STEP)
  960 FORMAT(14H DGSTPD       ,5(E20.13,3X))
      END
      SUBROUTINE DJAIRY(X,RX,C,AI,DAI)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DJAIRY
C***REFER TO  DBESJ,DBESY
C
C     1-2-74
C                  DJAIRY computes the Airy function AI(X)
C                   and its derivative DAI(X) for DASYJY
C
C                                   INPUT
C
C         X - Argument, computed by DASYJY, X unrestricted
C        RX - RX=DSQRT(DABS(X)), computed by DASYJY
C         C - C=2.*(DABS(X)**1.5)/3., computed by DASYJY
C
C                                  OUTPUT
C
C        AI - Value of function AI(X)
C       DAI - Value of the derivative DAI(X)
C
C                                Written by
C
C                                D. E. Amos
C                               S. L. Daniel
C                               M. K. Weston
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DJAIRY
C
      INTEGER I, J, M1, M1D, M2, M2D, M3, M3D, M4, M4D, N1, N1D, N2,
     1 N2D, N3, N3D, N4, N4D
      DOUBLE PRECISION A,AI,AJN,AJP,AK1,AK2,AK3,B,C,CCV,CON1,CON2,
     1 CON3, CON4, CON5, CV, DA, DAI, DAJN, DAJP, DAK1, DAK2, DAK3,
     2 DB, EC, E1, E2, FPI12, F1, F2, RTRX, RX, SCV, T, TEMP1, TEMP2,
     3 TT, X
      DIMENSION AJP(19), AJN(19), A(15), B(15)
      DIMENSION AK1(14), AK2(23), AK3(14)
      DIMENSION DAJP(19), DAJN(19), DA(15), DB(15)
      DIMENSION DAK1(14), DAK2(24), DAK3(14)
      SAVE N1, N2, N3, N4, M1, M2, M3, M4, FPI12, CON1, CON2, CON3,
     1 CON4, CON5, AK1, AK2, AK3, AJP, AJN, A, B,
     2 N1D, N2D, N3D, N4D, M1D, M2D, M3D, M4D, DAK1, DAK2, DAK3,
     3 DAJP, DAJN, DA, DB
      DATA N1,N2,N3,N4/14,23,19,15/
      DATA M1,M2,M3,M4/12,21,17,13/
      DATA FPI12,CON1,CON2,CON3,CON4,CON5/
     1 1.30899693899575D+00, 6.66666666666667D-01, 5.03154716196777D+00,
     2 3.80004589867293D-01, 8.33333333333333D-01, 8.66025403784439D-01/
      DATA AK1(1), AK1(2), AK1(3), AK1(4), AK1(5), AK1(6), AK1(7),
     1     AK1(8), AK1(9), AK1(10),AK1(11),AK1(12),AK1(13),
     2     AK1(14)         / 2.20423090987793D-01,-1.25290242787700D-01,
     3 1.03881163359194D-02, 8.22844152006343D-04,-2.34614345891226D-04,
     4 1.63824280172116D-05, 3.06902589573189D-07,-1.29621999359332D-07,
     5 8.22908158823668D-09, 1.53963968623298D-11,-3.39165465615682D-11,
     6 2.03253257423626D-12,-1.10679546097884D-14,-5.16169497785080D-15/
      DATA AK2(1), AK2(2), AK2(3), AK2(4), AK2(5), AK2(6), AK2(7),
     1     AK2(8), AK2(9), AK2(10),AK2(11),AK2(12),AK2(13),AK2(14),
     2     AK2(15),AK2(16),AK2(17),AK2(18),AK2(19),AK2(20),AK2(21),
     3     AK2(22),AK2(23) / 2.74366150869598D-01, 5.39790969736903D-03,
     4-1.57339220621190D-03, 4.27427528248750D-04,-1.12124917399925D-04,
     5 2.88763171318904D-05,-7.36804225370554D-06, 1.87290209741024D-06,
     6-4.75892793962291D-07, 1.21130416955909D-07,-3.09245374270614D-08,
     7 7.92454705282654D-09,-2.03902447167914D-09, 5.26863056595742D-10,
     8-1.36704767639569D-10, 3.56141039013708D-11,-9.31388296548430D-12,
     9 2.44464450473635D-12,-6.43840261990955D-13, 1.70106030559349D-13,
     1-4.50760104503281D-14, 1.19774799164811D-14,-3.19077040865066D-15/
      DATA AK3(1), AK3(2), AK3(3), AK3(4), AK3(5), AK3(6), AK3(7),
     1     AK3(8), AK3(9), AK3(10),AK3(11),AK3(12),AK3(13),
     2     AK3(14)         / 2.80271447340791D-01,-1.78127042844379D-03,
     3 4.03422579628999D-05,-1.63249965269003D-06, 9.21181482476768D-08,
     4-6.52294330229155D-09, 5.47138404576546D-10,-5.24408251800260D-11,
     5 5.60477904117209D-12,-6.56375244639313D-13, 8.31285761966247D-14,
     6-1.12705134691063D-14, 1.62267976598129D-15,-2.46480324312426D-16/
      DATA AJP(1), AJP(2), AJP(3), AJP(4), AJP(5), AJP(6), AJP(7),
     1     AJP(8), AJP(9), AJP(10),AJP(11),AJP(12),AJP(13),AJP(14),
     2     AJP(15),AJP(16),AJP(17),AJP(18),
     3     AJP(19)         / 7.78952966437581D-02,-1.84356363456801D-01,
     4 3.01412605216174D-02, 3.05342724277608D-02,-4.95424702513079D-03,
     5-1.72749552563952D-03, 2.43137637839190D-04, 5.04564777517082D-05,
     6-6.16316582695208D-06,-9.03986745510768D-07, 9.70243778355884D-08,
     7 1.09639453305205D-08,-1.04716330588766D-09,-9.60359441344646D-11,
     8 8.25358789454134D-12, 6.36123439018768D-13,-4.96629614116015D-14,
     9-3.29810288929615D-15, 2.35798252031104D-16/
      DATA AJN(1), AJN(2), AJN(3), AJN(4), AJN(5), AJN(6), AJN(7),
     1     AJN(8), AJN(9), AJN(10),AJN(11),AJN(12),AJN(13),AJN(14),
     2     AJN(15),AJN(16),AJN(17),AJN(18),
     3     AJN(19)         / 3.80497887617242D-02,-2.45319541845546D-01,
     4 1.65820623702696D-01, 7.49330045818789D-02,-2.63476288106641D-02,
     5-5.92535597304981D-03, 1.44744409589804D-03, 2.18311831322215D-04,
     6-4.10662077680304D-05,-4.66874994171766D-06, 7.15218807277160D-07,
     7 6.52964770854633D-08,-8.44284027565946D-09,-6.44186158976978D-10,
     8 7.20802286505285D-11, 4.72465431717846D-12,-4.66022632547045D-13,
     9-2.67762710389189D-14, 2.36161316570019D-15/
      DATA A(1),   A(2),   A(3),   A(4),   A(5),   A(6),   A(7),
     1     A(8),   A(9),   A(10),  A(11),  A(12),  A(13),  A(14),
     2     A(15)           / 4.90275424742791D-01, 1.57647277946204D-03,
     3-9.66195963140306D-05, 1.35916080268815D-07, 2.98157342654859D-07,
     4-1.86824767559979D-08,-1.03685737667141D-09, 3.28660818434328D-10,
     5-2.57091410632780D-11,-2.32357655300677D-12, 9.57523279048255D-13,
     6-1.20340828049719D-13,-2.90907716770715D-15, 4.55656454580149D-15,
     7-9.99003874810259D-16/
      DATA B(1),   B(2),   B(3),   B(4),   B(5),   B(6),   B(7),
     1     B(8),   B(9),   B(10),  B(11),  B(12),  B(13),  B(14),
     2     B(15)           / 2.78593552803079D-01,-3.52915691882584D-03,
     3-2.31149677384994D-05, 4.71317842263560D-06,-1.12415907931333D-07,
     4-2.00100301184339D-08, 2.60948075302193D-09,-3.55098136101216D-11,
     5-3.50849978423875D-11, 5.83007187954202D-12,-2.04644828753326D-13,
     6-1.10529179476742D-13, 2.87724778038775D-14,-2.88205111009939D-15,
     7-3.32656311696166D-16/
      DATA N1D,N2D,N3D,N4D/14,24,19,15/
      DATA M1D,M2D,M3D,M4D/12,22,17,13/
      DATA DAK1(1), DAK1(2), DAK1(3), DAK1(4), DAK1(5), DAK1(6),
     1     DAK1(7), DAK1(8), DAK1(9), DAK1(10),DAK1(11),DAK1(12),
     2    DAK1(13),DAK1(14)/ 2.04567842307887D-01,-6.61322739905664D-02,
     3-8.49845800989287D-03, 3.12183491556289D-03,-2.70016489829432D-04,
     4-6.35636298679387D-06, 3.02397712409509D-06,-2.18311195330088D-07,
     5-5.36194289332826D-10, 1.13098035622310D-09,-7.43023834629073D-11,
     6 4.28804170826891D-13, 2.23810925754539D-13,-1.39140135641182D-14/
      DATA DAK2(1), DAK2(2), DAK2(3), DAK2(4), DAK2(5), DAK2(6),
     1     DAK2(7), DAK2(8), DAK2(9), DAK2(10),DAK2(11),DAK2(12),
     2     DAK2(13),DAK2(14),DAK2(15),DAK2(16),DAK2(17),DAK2(18),
     3     DAK2(19),DAK2(20),DAK2(21),DAK2(22),DAK2(23),
     4     DAK2(24)        / 2.93332343883230D-01,-8.06196784743112D-03,
     5 2.42540172333140D-03,-6.82297548850235D-04, 1.85786427751181D-04,
     6-4.97457447684059D-05, 1.32090681239497D-05,-3.49528240444943D-06,
     7 9.24362451078835D-07,-2.44732671521867D-07, 6.49307837648910D-08,
     8-1.72717621501538D-08, 4.60725763604656D-09,-1.23249055291550D-09,
     9 3.30620409488102D-10,-8.89252099772401D-11, 2.39773319878298D-11,
     1-6.48013921153450D-12, 1.75510132023731D-12,-4.76303829833637D-13,
     2 1.29498241100810D-13,-3.52679622210430D-14, 9.62005151585923D-15,
     3-2.62786914342292D-15/
      DATA DAK3(1), DAK3(2), DAK3(3), DAK3(4), DAK3(5), DAK3(6),
     1     DAK3(7), DAK3(8), DAK3(9), DAK3(10),DAK3(11),DAK3(12),
     2    DAK3(13),DAK3(14)/ 2.84675828811349D-01, 2.53073072619080D-03,
     3-4.83481130337976D-05, 1.84907283946343D-06,-1.01418491178576D-07,
     4 7.05925634457153D-09,-5.85325291400382D-10, 5.56357688831339D-11,
     5-5.90889094779500D-12, 6.88574353784436D-13,-8.68588256452194D-14,
     6 1.17374762617213D-14,-1.68523146510923D-15, 2.55374773097056D-16/
      DATA DAJP(1), DAJP(2), DAJP(3), DAJP(4), DAJP(5), DAJP(6),
     1     DAJP(7), DAJP(8), DAJP(9), DAJP(10),DAJP(11),DAJP(12),
     2     DAJP(13),DAJP(14),DAJP(15),DAJP(16),DAJP(17),DAJP(18),
     3     DAJP(19)        / 6.53219131311457D-02,-1.20262933688823D-01,
     4 9.78010236263823D-03, 1.67948429230505D-02,-1.97146140182132D-03,
     5-8.45560295098867D-04, 9.42889620701976D-05, 2.25827860945475D-05,
     6-2.29067870915987D-06,-3.76343991136919D-07, 3.45663933559565D-08,
     7 4.29611332003007D-09,-3.58673691214989D-10,-3.57245881361895D-11,
     8 2.72696091066336D-12, 2.26120653095771D-13,-1.58763205238303D-14,
     9-1.12604374485125D-15, 7.31327529515367D-17/
      DATA DAJN(1), DAJN(2), DAJN(3), DAJN(4), DAJN(5), DAJN(6),
     1     DAJN(7), DAJN(8), DAJN(9), DAJN(10),DAJN(11),DAJN(12),
     2     DAJN(13),DAJN(14),DAJN(15),DAJN(16),DAJN(17),DAJN(18),
     3     DAJN(19)        / 1.08594539632967D-02, 8.53313194857091D-02,
     4-3.15277068113058D-01,-8.78420725294257D-02, 5.53251906976048D-02,
     5 9.41674060503241D-03,-3.32187026018996D-03,-4.11157343156826D-04,
     6 1.01297326891346D-04, 9.87633682208396D-06,-1.87312969812393D-06,
     7-1.50798500131468D-07, 2.32687669525394D-08, 1.59599917419225D-09,
     8-2.07665922668385D-10,-1.24103350500302D-11, 1.39631765331043D-12,
     9 7.39400971155740D-14,-7.32887475627500D-15/
      DATA DA(1),  DA(2),  DA(3),  DA(4),  DA(5),  DA(6),  DA(7),
     1     DA(8),  DA(9),  DA(10), DA(11), DA(12), DA(13), DA(14),
     2     DA(15)          / 4.91627321104601D-01, 3.11164930427489D-03,
     3 8.23140762854081D-05,-4.61769776172142D-06,-6.13158880534626D-08,
     4 2.87295804656520D-08,-1.81959715372117D-09,-1.44752826642035D-10,
     5 4.53724043420422D-11,-3.99655065847223D-12,-3.24089119830323D-13,
     6 1.62098952568741D-13,-2.40765247974057D-14, 1.69384811284491D-16,
     7 8.17900786477396D-16/
      DATA DB(1),  DB(2),  DB(3),  DB(4),  DB(5),  DB(6),  DB(7),
     1     DB(8),  DB(9),  DB(10), DB(11), DB(12), DB(13), DB(14),
     2     DB(15)          /-2.77571356944231D-01, 4.44212833419920D-03,
     3-8.42328522190089D-05,-2.58040318418710D-06, 3.42389720217621D-07,
     4-6.24286894709776D-09,-2.36377836844577D-09, 3.16991042656673D-10,
     5-4.40995691658191D-12,-5.18674221093575D-12, 9.64874015137022D-13,
     6-4.90190576608710D-14,-1.77253430678112D-14, 5.55950610442662D-15,
     7-7.11793337579530D-16/
C***FIRST EXECUTABLE STATEMENT  DJAIRY
      IF (X.LT.0.0D0) GO TO 90
      IF (C.GT.5.0D0) GO TO 60
      IF (X.GT.1.20D0) GO TO 30
      T = (X+X-1.2D0)*CON4
      TT = T + T
      J = N1
      F1 = AK1(J)
      F2 = 0.0D0
      DO 10 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK1(J)
        F2 = TEMP1
   10 CONTINUE
      AI = T*F1 - F2 + AK1(1)
C
      J = N1D
      F1 = DAK1(J)
      F2 = 0.0D0
      DO 20 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK1(J)
        F2 = TEMP1
   20 CONTINUE
      DAI = -(T*F1-F2+DAK1(1))
      RETURN
C
   30 CONTINUE
      T = (X+X-CON2)*CON3
      TT = T + T
      J = N2
      F1 = AK2(J)
      F2 = 0.0D0
      DO 40 I=1,M2
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK2(J)
        F2 = TEMP1
   40 CONTINUE
      RTRX = DSQRT(RX)
      EC = DEXP(-C)
      AI = EC*(T*F1-F2+AK2(1))/RTRX
      J = N2D
      F1 = DAK2(J)
      F2 = 0.0D0
      DO 50 I=1,M2D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK2(J)
        F2 = TEMP1
   50 CONTINUE
      DAI = -EC*(T*F1-F2+DAK2(1))*RTRX
      RETURN
C
   60 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N1
      F1 = AK3(J)
      F2 = 0.0D0
      DO 70 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK3(J)
        F2 = TEMP1
   70 CONTINUE
      RTRX = DSQRT(RX)
      EC = DEXP(-C)
      AI = EC*(T*F1-F2+AK3(1))/RTRX
      J = N1D
      F1 = DAK3(J)
      F2 = 0.0D0
      DO 80 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK3(J)
        F2 = TEMP1
   80 CONTINUE
      DAI = -RTRX*EC*(T*F1-F2+DAK3(1))
      RETURN
C
   90 CONTINUE
      IF (C.GT.5.0D0) GO TO 120
      T = 0.4D0*C - 1.0D0
      TT = T + T
      J = N3
      F1 = AJP(J)
      E1 = AJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 100 I=1,M3
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + AJP(J)
        E1 = TT*E1 - E2 + AJN(J)
        F2 = TEMP1
        E2 = TEMP2
  100 CONTINUE
      AI = (T*E1-E2+AJN(1)) - X*(T*F1-F2+AJP(1))
      J = N3D
      F1 = DAJP(J)
      E1 = DAJN(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 110 I=1,M3D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DAJP(J)
        E1 = TT*E1 - E2 + DAJN(J)
        F2 = TEMP1
        E2 = TEMP2
  110 CONTINUE
      DAI = X*X*(T*F1-F2+DAJP(1)) + (T*E1-E2+DAJN(1))
      RETURN
C
  120 CONTINUE
      T = 10.0D0/C - 1.0D0
      TT = T + T
      J = N4
      F1 = A(J)
      E1 = B(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 130 I=1,M4
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + A(J)
        E1 = TT*E1 - E2 + B(J)
        F2 = TEMP1
        E2 = TEMP2
  130 CONTINUE
      TEMP1 = T*F1 - F2 + A(1)
      TEMP2 = T*E1 - E2 + B(1)
      RTRX = DSQRT(RX)
      CV = C - FPI12
      CCV = DCOS(CV)
      SCV = DSIN(CV)
      AI = (TEMP1*CCV-TEMP2*SCV)/RTRX
      J = N4D
      F1 = DA(J)
      E1 = DB(J)
      F2 = 0.0D0
      E2 = 0.0D0
      DO 140 I=1,M4D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DA(J)
        E1 = TT*E1 - E2 + DB(J)
        F2 = TEMP1
        E2 = TEMP2
  140 CONTINUE
      TEMP1 = T*F1 - F2 + DA(1)
      TEMP2 = T*E1 - E2 + DB(1)
      E1 = CCV*CON5 + 0.5D0*SCV
      E2 = SCV*CON5 - 0.5D0*CCV
      DAI = (TEMP1*E1-TEMP2*E2)*RTRX
      RETURN
      END
      DOUBLE PRECISION FUNCTION DLNGAM(X)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DLNGAM
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C7A
C***KEYWORDS  GAMMA FUNCTION,LOGARITHM,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Computes the log of the absolute value of the Gamma
C            function
C***DESCRIPTION
C
C DLNGAM(X) computes the logarithm of the absolute value of the
C gamma function at X.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH,D9LGMC,DGAMMA,XERROR
C***END PROLOGUE  DLNGAM
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX,
     1  Y, DINT, DGAMMA, D9LGMC, D1MACH, TEMP
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA XMAX, DXREL / 2*0.D0 /
C***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (XMAX.NE.0.D0) GO TO 10
      TEMP = 1.D0/DLOG(D1MACH(2))
      XMAX = TEMP*D1MACH(2)
      DXREL = DSQRT (D1MACH(4))
C
 10   Y = DABS (X)
      IF (Y.GT.10.D0) GO TO 20
C
C DLOG (DABS (DGAMMA(X)) ) FOR DABS(X) .LE. 10.0
C
      DLNGAM = DLOG (DABS (DGAMMA(X)) )
      RETURN
C
C DLOG ( DABS (DGAMMA(X)) ) FOR DABS(X) .GT. 10.0
C
 20   IF (Y.GT.XMAX) CALL XERROR ( 'DLNGAM  DABS(X) SO BIG DLNGAM OVERFL
     1OWS', 39, 2, 2)
C
      IF (X.GT.0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*DLOG(X) - X + D9LGMC(Y)
      IF (X.GT.0.D0) RETURN
C
      SINPIY = DABS (DSIN(PI*Y))
      IF (SINPIY.EQ.0.D0) CALL XERROR ( 'DLNGAM  X IS A NEGATIVE INTEGER
     1', 31, 3, 2)
C
      IF (DABS ((X-DINT(X-0.5D0))/X).LT.DXREL) CALL XERROR ( 'DLNGAM  AN
     1SWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 68, 1
     2, 1)
C
      DLNGAM = SQPI2L + (X-0.5D0)*DLOG(Y) - X - DLOG(SINPIY) - D9LGMC(Y)
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION DNRM2(N,DX,INCX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DNRM2
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A3B
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SNRM2-S DNRM2-D SCNRM2-C),
C             EUCLIDEAN LENGTH,EUCLIDEAN NORM,L2,LINEAR ALGEBRA,UNITARY,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Euclidean length (L2 norm) of d.p. vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX() with storage
C     increment INCX .
C     If    N .LE. 0 return with result = 0.
C     If N .GE. 1 then INCX must be .GE. 1
C
C           C.L. Lawson, 1978 Jan 08
C
C     Four phase method     using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  DSQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  DSQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm..
C
C     Phase 1    scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI..
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as followS..
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
      SAVE CUTLO, CUTHI, ZERO, ONE
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DNRM2
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C***FIRST EXECUTABLE STATEMENT  DNRM2
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,
     1   MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,R,LR,QTF,WA1,
     2   WA2,WA3,WA4)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DNSQ
C***DATE WRITTEN   800301   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  F2A
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(SNSQ-S DNSQ-D),
C             NONLINEAR SQUARE SYSTEM,POWELL HYBRID METHOD,ZEROES,ZEROS
C***AUTHOR  HIEBERT K.L. (SNLA)
C***PURPOSE  To find a zero of a system of N nonlinear functions in
C            N variables by a modification of the Powell hybrid method.
C            This code is the combination of the MINPACK codes (Argonne)
C            HYBRD and HYBRDJ.
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNSQ is to find a zero of a system of N nonlinear
C       functions in N variables by a modification of the Powell
C       hybrid method.  The user must provide a subroutine which
C       calculates the functions.  The user has the option of either to
C       provide a subroutine which calculates the Jacobian or to let the
C       code calculate it by a forward-difference approximation.
C       This code is the combination of the MINPACK codes (Argonne)
C       HYBRD and HYBRDJ.
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
C      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
C      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
C       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR
C       DOUBLE PRECISION
C       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
C      *     WA1(N),WA2(N),WA3(N),WA4(N)
C       EXTERNAL FCN,JAC
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNSQ and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNSQ.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an EXTERNAL statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         EXTERNAL statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=1, then the user must supply the
C         Jacobian through the subroutine JAC.  If IOPT=2, then the
C         code will approximate the Jacobian by forward-differencing.
C
C       N is a positive integer input variable set to the number of
C         functions and variables.
C
C       X is an array of length N.  On input X must contain an initial
C         estimate of the solution vector.  On output X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       FJAC is an output N by N array which contains the orthogonal
C         matrix Q produced by the QR factorization of the final
C         approximate Jacobian.
C
C       LDFJAC is a positive integer input variable not less than N
C         which specifies the leading dimension of the array FJAC.
C
C       XTOL is a nonnegative input variable.  Termination occurs when
C         the relative error between two consecutive iterates is at most
C         XTOL.  Therefore, XTOL measures the relative error desired in
C         the approximate solution.  Section 4 contains more details
C         about XTOL.
C
C       MAXFEV is a positive integer input variable.  Termination occurs
C         when the number of calls to FCN is at least MAXFEV by the end
C         of an iteration.
C
C       ML is a nonnegative integer input variable which specifies the
C         number of subdiagonals within the band of the Jacobian matrix.
C         If the Jacobian is not banded or IOPT=1, set ML to at
C         least N - 1.
C
C       MU is a nonnegative integer input variable which specifies the
C         number of superdiagonals within the band of the Jacobian
C         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
C         least N - 1.
C
C       EPSFCN is an input variable used in determining a suitable step
C         for the forward-difference approximation.  This approximation
C         assumes that the relative errors in the functions are of the
C         order of EPSFCN.  If EPSFCN is less than the machine
C         precision, it is assumed that the relative errors in the
C         functions are of the order of the machine precision.  If
C         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
C         argument).
C
C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
C         internally set.  If MODE = 2, DIAG must contain positive
C         entries that serve as implicit (multiplicative) scale factors
C         for the variables.
C
C       MODE is an integer input variable.  If MODE = 1, the variables
C         will be scaled internally.  If MODE = 2, the scaling is
C         specified by the input DIAG.  Other values of MODE are
C         equivalent to MODE = 1.
C
C       FACTOR is a positive input variable used in determining the
C         initial step bound.  This bound is set to the product of
C         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
C         FACTOR itself.  In most cases FACTOR should lie in the
C         interval (.1,100.).  100. is a generally recommended value.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. appropriate
C         print statements must be added to FCN(see example).  If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  Improper input parameters.
C
C         INFO = 1  Relative error between two consecutive iterates is
C                   at most XTOL.
C
C         INFO = 2  Number of calls to FCN has reached or exceeded
C                   MAXFEV.
C
C         INFO = 3  XTOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  Iteration is not making good progress, as measured
C                   by the improvement from the last five Jacobian
C                   evaluations.
C
C         INFO = 5  Iteration is not making good progress, as measured
C                   by the improvement from the last ten iterations.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       NFEV is an integer output variable set to the number of calls to
C         FCN.
C
C       NJEV is an integer output variable set to the number of calls to
C         JAC. (If IOPT=2, then NJEV is set to zero.)
C
C       R is an output array of length LR which contains the upper
C         triangular matrix produced by the QR factorization of the
C         final approximate Jacobian, stored rowwise.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       QTF is an output array of length N which contains the vector
C         (Q transpose)*FVEC.
C
C       WA1, WA2, WA3, and WA4 are work arrays of length N.
C
C
C 4. Successful completion.
C
C       The accuracy of DNSQ is controlled by the convergence parameter
C       XTOL.  This parameter is used in a test which makes a comparison
C       between the approximation X and a solution XSOL.  DNSQ
C       terminates when the test is satisfied.  If the convergence
C       parameter is less than the machine precision (as defined by the
C       function D1MACH(4)), then DNSQ only attempts to satisfy the test
C       defined by the machine precision.  Further progress is not
C       usually possible.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then DNSQ may incorrectly indicate
C       convergence.  The coding of the Jacobian can be checked by the
C       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
C       then the validity of the answer can be checked, for example, by
C       rerunning DNSQ with a tighter tolerance.
C
C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
C         vector Z and D is the diagonal matrix whose entries are
C         defined by the array DIAG, then this test attempts to
C         guarantee that
C
C               DENORM(D*(X-XSOL)) .LE. XTOL*DENORM(D*XSOL).
C
C         If this condition is satisfied with XTOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller
C         components of D*X may have large relative errors, but the fast
C         rate of convergence of DNSQ usually avoids this possibility.
C         Unless high precision solutions are required, the recommended
C         value for XTOL is the square root of the machine precision.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNSQ can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of
C       function evaluations, or lack of good progress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
C         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or
C         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0,
C         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNSQ.  In this
C         case, it may be possible to remedy the situation by rerunning
C         DNSQ with a smaller value of FACTOR.
C
C       Excessive Number of Function Evaluations.  A reasonable value
C         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
C         If the number of calls to FCN reaches MAXFEV, then this
C         indicates that the routine is converging very slowly as
C         measured by the progress of FVEC, and INFO is set to 2. This
C         situation should be unusual because, as indicated below, lack
C         of good progress is usually diagnosed earlier by DNSQ,
C         causing termination with info = 4 or INFO = 5.
C
C       Lack of Good Progress.  DNSQ searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this
C         situation, the iteration eventually fails to make good
C         progress.  In particular, this will happen if the system does
C         not have a zero.  If the system has a zero, rerunning DNSQ
C         from a different starting point may be helpful.
C
C
C 6. Characteristics of The Algorithm.
C
C       DNSQ is a modification of the Powell Hybrid method.  Two of its
C       main characteristics involve the choice of the correction as a
C       convex combination of the Newton and scaled gradient directions,
C       and the updating of the Jacobian by the rank-1 method of
C       Broyden.  The choice of the correction guarantees (under
C       reasonable conditions) global convergence for starting points
C       far from the solution and a fast rate of convergence.  The
C       Jacobian is calculated at the starting point by either the
C       user-supplied subroutine or a forward-difference approximation,
C       but it is not recalculated until the rank-1 method fails to
C       produce satisfactory progress.
C
C       Timing.  The time required by DNSQ to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by DNSQ is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of DNSQ will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C   **For a completely set up example problem see the LONG DESCRIPTION**
C***LONG DESCRIPTION
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C C     **********
C
C       PROGRAM TEST(INPUT,OUTPUT,TAPE6=OUTPUT)
C C
C C     Driver for DNSQ example.
C C
C       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
C      *        NWRITE
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM
C       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
C      *     WA1(9),WA2(9),WA3(9),WA4(9)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C C
C       LDFJAC = 9
C       LR = 45
C C
C C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C C     THIS IS THE RECOMMENDED SETTING.
C C
C       XTOL = DSQRT(D1MACH(4))
C C
C       MAXFEV = 2000
C       ML = 1
C       MU = 1
C       EPSFCN = 0.E0
C       MODE = 2
C       DO 20 J = 1, 9
C          DIAG(J) = 1.E0
C    20    CONTINUE
C       FACTOR = 1.E2
C       NPRINT = 0
C C
C       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
C      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
C      *           R,LR,QTF,WA1,WA2,WA3,WA4)
C       FNORM = DENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(N)
C       INTEGER K
C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 K = 1, N
C          TEMP = (THREE - TWO*X(K))*X(K)
C          TEMP1 = ZERO
C          IF (K .NE. 1) TEMP1 = X(K-1)
C          TEMP2 = ZERO
C          IF (K .NE. N) TEMP2 = X(K+1)
C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
C    10    CONTINUE
C       RETURN
C       END
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       Final L2 norm of the residuals  0.1192636E-07
C
C       Number of function evaluations        14
C
C       Exit parameter                         1
C
C       Final approximate solution
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C***REFERENCES  POWELL, M. J. D., *A HYBRID METHOD FOR NONLINEAR
C                 EQUATIONS*, 'NUMERICAL METHODS FOR NONLINEAR
C                 ALGEBRAIC EQUATIONS', P. RABINOWITZ, EDITOR,
C                 GORDON AND BREACH, 1970.
C***ROUTINES CALLED  D1MACH,D1MPYQ,D1UPDT,DDOGLG,DENORM,DFDJC1,DQFORM,
C                    DQRFAC,XERROR
C***END PROLOGUE  DNSQ
      INTEGER MIN0,MOD
      DOUBLE PRECISION D1MACH,DENORM
      DOUBLE PRECISION DMAX1,DMIN1,DABS
      INTEGER I, IFLAG, INFO, IOPT, ITER, IWA(1), J, JM1, L, LDFJAC,
     1     LR, MAXFEV, ML, MODE, MU, N, NCFAIL, NCSUC, NFEV, NJEV,
     2     NPRINT, NSLOW1, NSLOW2
      DOUBLE PRECISION ACTRED, DELTA, DIAG(N), EPSFCN, EPSMCH, FACTOR,
     1     FJAC(LDFJAC,N), FNORM, FNORM1, FVEC(N), ONE, P0001, P001,
     2     P1, P5, PNORM, PRERED, QTF(N), R(LR), RATIO, SUM, TEMP,
     3     WA1(N), WA2(N), WA3(N), WA4(N), X(N), XNORM, XTOL, ZERO
      EXTERNAL FCN
      LOGICAL JEVAL,SING
      SAVE ONE, P1, P5, P001, P0001, ZERO
      DATA ONE,P1,P5,P001,P0001,ZERO
     1     /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 320
C***FIRST EXECUTABLE STATEMENT  DNSQ
         EPSMCH = D1MACH(4)
C
         INFO = 0
         IFLAG = 0
         NFEV = 0
         NJEV = 0
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C
C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1       .OR. XTOL .LT. ZERO .OR. MAXFEV .LE. 0 .OR. ML .LT. 0
     2       .OR. MU .LT. 0 .OR. FACTOR .LE. ZERO .OR. LDFJAC .LT. N
     3       .OR. LR .LT. (N*(N + 1))/2) GO TO 320
         IF (MODE .NE. 2) GO TO 20
            DO 10 J = 1, N
C     .........EXIT
               IF (DIAG(J) .LE. ZERO) GO TO 320
   10       CONTINUE
   20    CONTINUE
C
C        EVALUATE THE FUNCTION AT THE STARTING POINT
C        AND CALCULATE ITS NORM.
C
         IFLAG = 1
         CALL FCN(N,X,FVEC,IFLAG)
         NFEV = 1
C     ...EXIT
         IF (IFLAG .LT. 0) GO TO 320
         FNORM = DENORM(N,FVEC)
C
C        INITIALIZE ITERATION COUNTER AND MONITORS.
C
         ITER = 1
         NCSUC = 0
         NCFAIL = 0
         NSLOW1 = 0
         NSLOW2 = 0
C
C        BEGINNING OF THE OUTER LOOP.
C
   30    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 90
               JEVAL = .TRUE.
C
C              CALCULATE THE JACOBIAN MATRIX.
C
               IF (IOPT .EQ. 2) GO TO 40
C
C                 USER SUPPLIES JACOBIAN
C
                  CALL JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
                  NJEV = NJEV + 1
               GO TO 50
   40          CONTINUE
C
C                 CODE APPROXIMATES THE JACOBIAN
C
                  IFLAG = 2
                  CALL DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,
     1                        EPSFCN,WA1,WA2)
                  NFEV = NFEV + MIN0(ML+MU+1,N)
   50          CONTINUE
C
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
C
C              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
               CALL DQRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C
C              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
C           ...EXIT
               IF (ITER .NE. 1) GO TO 90
               IF (MODE .EQ. 2) GO TO 70
                  DO 60 J = 1, N
                     DIAG(J) = WA2(J)
                     IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   60             CONTINUE
   70          CONTINUE
C
C              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
C              X AND INITIALIZE THE STEP BOUND DELTA.
C
               DO 80 J = 1, N
                  WA3(J) = DIAG(J)*X(J)
   80          CONTINUE
               XNORM = DENORM(N,WA3)
               DELTA = FACTOR*XNORM
               IF (DELTA .EQ. ZERO) DELTA = FACTOR
   90       CONTINUE
C
C           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
C
            DO 100 I = 1, N
               QTF(I) = FVEC(I)
  100       CONTINUE
            DO 140 J = 1, N
               IF (FJAC(J,J) .EQ. ZERO) GO TO 130
                  SUM = ZERO
                  DO 110 I = J, N
                     SUM = SUM + FJAC(I,J)*QTF(I)
  110             CONTINUE
                  TEMP = -SUM/FJAC(J,J)
                  DO 120 I = J, N
                     QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
C
C           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
C
            SING = .FALSE.
            DO 170 J = 1, N
               L = J
               JM1 = J - 1
               IF (JM1 .LT. 1) GO TO 160
               DO 150 I = 1, JM1
                  R(L) = FJAC(I,J)
                  L = L + N - I
  150          CONTINUE
  160          CONTINUE
               R(L) = WA1(J)
               IF (WA1(J) .EQ. ZERO) SING = .TRUE.
  170       CONTINUE
C
C           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
C
            CALL DQFORM(N,N,FJAC,LDFJAC,WA1)
C
C           RESCALE IF NECESSARY.
C
            IF (MODE .EQ. 2) GO TO 190
               DO 180 J = 1, N
                  DIAG(J) = DMAX1(DIAG(J),WA2(J))
  180          CONTINUE
  190       CONTINUE
C
C           BEGINNING OF THE INNER LOOP.
C
  200       CONTINUE
C
C              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
               IF (NPRINT .LE. 0) GO TO 210
                  IFLAG = 0
                  IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1               CALL FCN(N,X,FVEC,IFLAG)
C     ............EXIT
                  IF (IFLAG .LT. 0) GO TO 320
  210          CONTINUE
C
C              DETERMINE THE DIRECTION P.
C
               CALL DDOGLG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
               DO 220 J = 1, N
                  WA1(J) = -WA1(J)
                  WA2(J) = X(J) + WA1(J)
                  WA3(J) = DIAG(J)*WA1(J)
  220          CONTINUE
               PNORM = DENORM(N,WA3)
C
C              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
               IF (ITER .EQ. 1) DELTA = DMIN1(DELTA,PNORM)
C
C              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
               IFLAG = 1
               CALL FCN(N,WA2,WA4,IFLAG)
               NFEV = NFEV + 1
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
               FNORM1 = DENORM(N,WA4)
C
C              COMPUTE THE SCALED ACTUAL REDUCTION.
C
               ACTRED = -ONE
               IF (FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C              COMPUTE THE SCALED PREDICTED REDUCTION.
C
               L = 1
               DO 240 I = 1, N
                  SUM = ZERO
                  DO 230 J = I, N
                     SUM = SUM + R(L)*WA1(J)
                     L = L + 1
  230             CONTINUE
                  WA3(I) = QTF(I) + SUM
  240          CONTINUE
               TEMP = DENORM(N,WA3)
               PRERED = ZERO
               IF (TEMP .LT. FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C              REDUCTION.
C
               RATIO = ZERO
               IF (PRERED .GT. ZERO) RATIO = ACTRED/PRERED
C
C              UPDATE THE STEP BOUND.
C
               IF (RATIO .GE. P1) GO TO 250
                  NCSUC = 0
                  NCFAIL = NCFAIL + 1
                  DELTA = P5*DELTA
               GO TO 260
  250          CONTINUE
                  NCFAIL = 0
                  NCSUC = NCSUC + 1
                  IF (RATIO .GE. P5 .OR. NCSUC .GT. 1)
     1               DELTA = DMAX1(DELTA,PNORM/P5)
                  IF (DABS(RATIO-ONE) .LE. P1) DELTA = PNORM/P5
  260          CONTINUE
C
C              TEST FOR SUCCESSFUL ITERATION.
C
               IF (RATIO .LT. P0001) GO TO 280
C
C                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
                  DO 270 J = 1, N
                     X(J) = WA2(J)
                     WA2(J) = DIAG(J)*X(J)
                     FVEC(J) = WA4(J)
  270             CONTINUE
                  XNORM = DENORM(N,WA2)
                  FNORM = FNORM1
                  ITER = ITER + 1
  280          CONTINUE
C
C              DETERMINE THE PROGRESS OF THE ITERATION.
C
               NSLOW1 = NSLOW1 + 1
               IF (ACTRED .GE. P001) NSLOW1 = 0
               IF (JEVAL) NSLOW2 = NSLOW2 + 1
               IF (ACTRED .GE. P1) NSLOW2 = 0
C
C              TEST FOR CONVERGENCE.
C
               IF (DELTA .LE. XTOL*XNORM .OR. FNORM .EQ. ZERO) INFO = 1
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
               IF (NFEV .GE. MAXFEV) INFO = 2
               IF (P1*DMAX1(P1*DELTA,PNORM) .LE. EPSMCH*XNORM) INFO = 3
               IF (NSLOW2 .EQ. 5) INFO = 4
               IF (NSLOW1 .EQ. 10) INFO = 5
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              CRITERION FOR RECALCULATING JACOBIAN
C
C           ...EXIT
               IF (NCFAIL .EQ. 2) GO TO 310
C
C              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
C              AND UPDATE QTF IF NECESSARY.
C
               DO 300 J = 1, N
                  SUM = ZERO
                  DO 290 I = 1, N
                     SUM = SUM + FJAC(I,J)*WA4(I)
  290             CONTINUE
                  WA2(J) = (SUM - WA3(J))/PNORM
                  WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
                  IF (RATIO .GE. P0001) QTF(J) = SUM
  300          CONTINUE
C
C              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
C
               CALL D1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)
               CALL D1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
               CALL D1MPYQ(1,N,QTF,1,WA2,WA3)
C
C              END OF THE INNER LOOP.
C
               JEVAL = .FALSE.
            GO TO 200
  310       CONTINUE
C
C           END OF THE OUTER LOOP.
C
         GO TO 30
  320 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(N,X,FVEC,IFLAG)
      IF (INFO .LT. 0)
     1   CALL XERROR( 'DNSQ   -- EXECUTION TERMINATED BECAUSE USER SET I
     2FLAG NEGATIVE.',63,1,1)
      IF (INFO .EQ. 0)
     1   CALL XERROR( 'DNSQ   -- INVALID INPUT PARAMETER.',34,2,1)
      IF (INFO .EQ. 2)
     1   CALL XERROR( 'DNSQ   -- TOO MANY FUNCTION EVALUATIONS.',40,9,1)
      IF (INFO .EQ. 3)
     1   CALL XERROR(  'DNSQ   -- XTOL TOO SMALL. NO FURTHER IMPROVEMENT
     2 POSSIBLE.',    58,3,1)
      IF (INFO .GT. 4)
     1   CALL XERROR( 'DNSQ   -- ITERATION NOT MAKING GOOD PROGRESS.',
     2               45,1,1)
      RETURN
      END
      SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DNSQE
C***DATE WRITTEN   800301   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  F2A
C***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(SNSQE-S DNSQE-D),
C             EASY-TO-USE,NONLINEAR SQUARE SYSTEM,POWELL HYBRID METHOD,
C             ZEROES,ZEROS
C***AUTHOR  HIEBERT K.L. (SNLA)
C***PURPOSE  The easy-to-use version of DNSQ which finds a zero of a
C            system of N nonlinear functions in N variables by a
C            modification of the Powell hybrid method.  This code is a
C            combination of the MINPACK codes HYBRD1 and HYBRJ1.
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNSQE is to find a zero of a system of N
C       nonlinear functions in N variables by a modification of the
C       Powell hybrid method.  This is done by using the more general
C       nonlinear equation solver DNSQ.  The user must provide a
C       subroutine which calculates the functions.  The user has the
C       option of either to provide a subroutine which calculates the
C       Jacobian or to let the code calculate it by a forward-difference
C       approximation.  This code is the combination of the MINPACK
C       codes (Argonne) HYBRD1 and HYBRJ1.
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
C      *                  WA,LWA)
C       INTEGER IOPT,N,NPRINT,INFO,LWA
C       DOUBLE PRECISION TOL
C       DOUBLE PRECISION X(N),FVEC(N),WA(LWA)
C       EXTERNAL FCN,JAC
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNSQE and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNSQE.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an external statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         Calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNSQE.  In this case set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         external statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of DNSQE. In this case set
C         IFLAG to a negative integer.
C
C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=1, then the user must supply the
C         Jacobian through the subroutine JAC.  If IOPT=2, then the
C         code will approximate the Jacobian by forward-differencing.
C
C       N is a positive integer input variable set to the number of
C         functions and variables.
C
C       X is an array of length N.  On input X must contain an initial
C         estimate of the solution vector.  On output X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       TOL is a nonnegative input variable.  Termination occurs when
C         the algorithm estimates that the relative error between X and
C         the solution is at most TOL.  Section 4 contains more details
C         about TOL.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN(see example).  If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  Improper input parameters.
C
C         INFO = 1  Algorithm estimates that the relative error between
C                   X and the solution is at most TOL.
C
C         INFO = 2  Number of calls to FCN has reached or exceeded
C                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
C
C         INFO = 3  TOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  Iteration is not making good progress.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       WA is a work array of length LWA.
C
C       LWA is a positive integer input variable not less than
C         (3*N**2+13*N))/2.
C
C 4. Successful Completion.
C
C       The accuracy of DNSQE is controlled by the convergence parameter
C       TOL.  This parameter is used in a test which makes a comparison
C       between the approximation X and a solution XSOL.  DNSQE
C       terminates when the test is satisfied.  If TOL is less than the
C       machine precision (as defined by the  function D1MACH(4)), then
C       DNSQE only attempts to satisfy the test defined by the machine
C       precision.  Further progress is not usually possible.  Unless
C       high precision solutions are required, the recommended value
C       for TOL is the square root of the machine precision.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently. If these conditions are
C       not satisfied, then DNSQE may incorrectly indicate convergence.
C       The coding of the Jacobian can be checked by the subroutine
C       DCKDER.  If the Jacobian is coded correctly or IOPT=2, then
C       the validity of the answer can be checked, for example, by
C       rerunning DNSQE with a tighter tolerance.
C
C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
C         vector Z, then this test attempts to guarantee that
C
C               DENORM(X-XSOL) .LE. TOL*DENORM(XSOL).
C
C         If this condition is satisfied with TOL = 10**(-K), then the
C         larger components of X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller
C         components of X may have large relative errors, but the fast
C         rate of convergence of DNSQE usually avoids this possibility.
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNSQE can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of
C       function evaluations, errors in the functions, or lack of good
C       progress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or
C         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or
C         LWA .LT. (3*N**2+13*N)/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNSQE.  In this
C         case, it may be possible to remedy the situation by not
C         evaluating the functions here, but instead setting the
C         components of FVEC to numbers that exceed those in the initial
C         FVEC.
C
C       Excessive Number of Function Evaluations.  If the number of
C         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
C         IOPT=2, then this indicates that the routine is converging
C         very slowly as measured by the progress of FVEC, and INFO is
C         set to 2.  This situation should be unusual because, as
C         indicated below, lack of good progress is usually diagnosed
C         earlier by DNSQE, causing termination with INFO = 4.
C
C       Errors In the Functions.  When IOPT=2, the choice of step length
C         in the forward-difference approximation to the Jacobian
C         assumes that the relative errors in the functions are of the
C         order of the machine precision.  If this is not the case,
C         DNSQE may fail (usually with INFO = 4).  The user should
C         then either use DNSQ and set the step length or use IOPT=1
C         and supply the Jacobian.
C
C       Lack of Good Progress.  DNSQE searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this
C         situation, the iteration eventually fails to make good
C         progress.  In particular, this will happen if the system does
C         not have a zero.  If the system has a zero, rerunning DNSQE
C         from a different starting point may be helpful.
C
C 6. Characteristics of The Algorithm.
C
C       DNSQE is a modification of the Powell Hybrid method.  Two of
C       its main characteristics involve the choice of the correction as
C       a convex combination of the Newton and scaled gradient
C       directions, and the updating of the Jacobian by the rank-1
C       method of Broyden.  The choice of the correction guarantees
C       (under reasonable conditions) global convergence for starting
C       points far from the solution and a fast rate of convergence.
C       The Jacobian is calculated at the starting point by either the
C       user-supplied subroutine or a forward-difference approximation,
C       but it is not recalculated until the rank-1 method fails to
C       produce satisfactory progress.
C
C       Timing.  The time required by DNSQE to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by DNSQE is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of DNSQE will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  DNSQE requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C  *For a completely set up example problem see the LONG DESCRIPTION**
C***LONG DESCRIPTION
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C
C       **********
C
C       PROGRAM TEST(INPUT,OUTPUT,TAPE6=OUTPUT)
C C
C C     DRIVER FOR DNSQE EXAMPLE.
C C
C       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
C       DOUBLE PRECISION TOL,FNORM
C       DOUBLE PRECISION X(9),FVEC(9),WA(180)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C
C       LWA = 180
C       NPRINT = 0
C C
C C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C C     THIS IS THE RECOMMENDED SETTING.
C C
C       TOL = DSQRT(D1MACH(4))
C C
C       CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
C       FNORM = DENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(N)
C       INTEGER K
C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
C       DO 10 K = 1, N
C          TEMP = (THREE - TWO*X(K))*X(K)
C          TEMP1 = ZERO
C          IF (K .NE. 1) TEMP1 = X(K-1)
C          TEMP2 = ZERO
C          IF (K .NE. N) TEMP2 = X(K+1)
C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
C    10    CONTINUE
C       RETURN
C       END
C
C       RESULTS OBTAINED WITH DIFFERENT COMPILERS OR MACHINES
C       MAY BE SLIGHTLY DIFFERENT.
C
C       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C***REFERENCES  POWELL, M. J. D., *A HYBRID METHOD FOR NONLINEAR
C                 EQUATIONS*, 'NUMERICAL METHODS FOR NONLINEAR
C                 ALGEBRAIC EQUATIONS', P. RABINOWITZ, EDITOR,
C                 GORDON AND BREACH, 1970.
C***ROUTINES CALLED  DNSQ,XERROR
C***END PROLOGUE  DNSQE
      INTEGER INDEX, INFO, IOPT, J, LR, LWA, MAXFEV, ML, MODE, MU, N,
     1     NFEV, NJEV, NPRINT
      DOUBLE PRECISION EPSFCN, FACTOR, FVEC(N), ONE, TOL, WA(LWA),
     1     X(N), XTOL, ZERO
      EXTERNAL FCN,JAC
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
C     BEGIN BLOCK PERMITTING ...EXITS TO 20
C***FIRST EXECUTABLE STATEMENT  DNSQE
         INFO = 0
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C
C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1       .OR. TOL .LT. ZERO .OR. LWA .LT. (3*N**2 + 13*N)/2)
     2      GO TO 20
C
C        CALL DNSQ.
C
         MAXFEV = 100*(N + 1)
         IF (IOPT .EQ. 2) MAXFEV = 2*MAXFEV
         XTOL = TOL
         ML = N - 1
         MU = N - 1
         EPSFCN = ZERO
         MODE = 2
         DO 10 J = 1, N
            WA(J) = ONE
   10    CONTINUE
         LR = (N*(N + 1))/2
         INDEX = 6*N + LR
         CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,WA(INDEX+1),N,XTOL,MAXFEV,ML,
     1             MU,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
     2             WA(6*N+1),LR,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),
     3             WA(5*N+1))
         IF (INFO .EQ. 5) INFO = 4
   20 CONTINUE
      IF (INFO .EQ. 0)
     1   CALL XERROR( 'DNSQE  -- INVALID INPUT PARAMETER.',34,2,1)
      RETURN
      END
      SUBROUTINE DPASB2(IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASB2
C***REFER TO  DCFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASB2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  DPASB2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB3(IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASB3
C***REFER TO  DCFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASB3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  DPASB3
      TAUR = -.5D0
      TAUI = .5D0*SQRT(3.0D0)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB4(IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASB4
C***REFER TO  DCFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASB4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  DPASB4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASB5
C***REFER TO  DCFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASB5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  DPASB5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(.1D0*PI)
      TI11 = SIN(.4D0*PI)
      TR12 = -SIN(.3D0*PI)
      TI12 = SIN(.2D0*PI)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF2(IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASF2
C***REFER TO  DCFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASF2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  DPASF2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
      DO 106 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF3(IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASF3
C***REFER TO  DCFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASF3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  DPASF3
      TAUR = -.5D0
      TAUI = -.5D0*SQRT(3.0D0)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF4(IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASF4
C***REFER TO  DCFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASF4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  DPASF4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASF5
C***REFER TO  DCFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASF5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  DPASF5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(.1D0*PI)
      TI11 = -SIN(.4D0*PI)
      TR12 = -SIN(.3D0*PI)
      TI12 = -SIN(.2D0*PI)
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 IF(IDO/2.LT.L1) GO TO 105
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
  105 DO 107 I=2,IDO,2
CDIR$ IVDEP
         DO 106 K=1,L1
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  106    CONTINUE
  107 CONTINUE
      RETURN
      END
      SUBROUTINE DPASSB(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASSB
C***REFER TO  DCFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASSB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
C***FIRST EXECUTABLE STATEMENT  DPASSB
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
CDIR$ IVDEP
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
CDIR$ IVDEP
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
CDIR$ IVDEP
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
CDIR$ IVDEP
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
CDIR$ IVDEP
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
CDIR$ IVDEP
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
CDIR$ IVDEP
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
CDIR$ IVDEP
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
CDIR$ IVDEP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
CDIR$ IVDEP
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
CDIR$ IVDEP
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
      SUBROUTINE DPASSF(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPASSF
C***REFER TO  DCFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DPASSF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
     2                CH2(IDL1,IP)
C***FIRST EXECUTABLE STATEMENT  DPASSF
      IDOT = IDO/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
CDIR$ IVDEP
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
CDIR$ IVDEP
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
CDIR$ IVDEP
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
CDIR$ IVDEP
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
CDIR$ IVDEP
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
CDIR$ IVDEP
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
CDIR$ IVDEP
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
CDIR$ IVDEP
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
CDIR$ IVDEP
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
CDIR$ IVDEP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
CDIR$ IVDEP
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
CDIR$ IVDEP
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DPCHDF(K,X,S,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHDF
C***REFER TO  DPCHCE,DPCHSP
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          DPCHDF:   DPCHIP Finite Difference Formula
C
C     Uses a divided difference formulation to compute a K-point approx-
C     imation to the derivative at X(K) based on the data in X and S.
C
C     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
C     derivative approximations.
C
C ----------------------------------------------------------------------
C
C     On input:
C        K      is the order of the desired derivative approximation.
C               K must be at least 3 (error return if not).
C        X      contains the K values of the independent variable.
C               X need not be ordered, but the values **MUST** be
C               distinct.  (Not checked here.)
C        S      contains the associated slope values:
C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
C               (Note that S need only be of length K-1.)
C
C     On return:
C        S      will be destroyed.
C        IERR   will be set to -1 if K.LT.2 .
C        DPCHDF  will be set to the desired derivative approximation if
C               IERR=0 or to zero if IERR=-1.
C
C ----------------------------------------------------------------------
C
C  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
C              Verlag (New York, 1978), pp. 10-16.
C
C***END PROLOGUE  DPCHDF
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHDF to PCHDF wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constant Zero to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER K, IERR
      DOUBLE PRECISION  X(K), S(K)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, J
      DOUBLE PRECISION  VALUE, ZERO
      DATA  ZERO /0.D0/
C
C  CHECK FOR LEGAL VALUE OF K.
C
C***FIRST EXECUTABLE STATEMENT  DPCHDF
      IF (K .LT. 3)  GO TO 5001
C
C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
C
C  EVALUATE DERIVATIVE AT X(K).
C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      DPCHDF = VALUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     K.LT.3 RETURN.
      IERR = -1
      CALL XERROR ('DPCHDF -- K LESS THAN THREE'
     *           , 0, IERR, 1)
      DPCHDF = ZERO
      RETURN
      END
      SUBROUTINE DPCHEV(N,X,F,D,NVAL,XVAL,FVAL,DVAL,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHEV
C***DATE WRITTEN   870828   (YYMMDD)
C***REVISION DATE  870828   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE
C             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the function and first derivative of a piecewise
C            cubic Hermite or spline function at an array of points
C            XVAL.  It is easy to use.
C***DESCRIPTION
C
C          DPCHEV:  Piecewise Cubic Hermite or Spline Derivative
C                   Evaluator. Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C                 Prentice Hall 1988
C
C     Evaluates the function and first derivative of the cubic Hermite
C     or spline function defined by  N, X, F, D, at the array of points
C     XVAL.
C
C
C     This is an easy to use driver for the routines by F.N. Fritsch
C     described in reference (2) below. Those also have other 
C     capabilities.
C
C ----------------------------------------------------------------------
C
C  Calling sequence: CALL  DPCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
C
C     INTEGER  N, NVAL, IERR
C     DOUBLE PRECISION  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) double precision array of independent variable 
C           values.  The elements of X must be strictly increasing:
C             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
C
C     F -- (input) double precision array of function values.  F(I) is
C           the value corresponding to X(I).
C
C     D -- (input) double precision array of derivative values.  
C          D(I) is the value corresponding to X(I).
C
C  NVAL -- (input) number of points at which the functions are to be
C           evaluated. ( Error return if NVAL.LT.1 )
C
C  XVAL -- (input) double precision array of points at which the 
C          functions are to be evaluated.
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XVAL are increasing relative to X;
C              that is,   XVAL(J) .GE. X(I)
C              implies    XVAL(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XVAL are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C  FVAL -- (output) double precision array of values of the cubic 
C          Hermite function defined by  N, X, F, D  at the points  XVAL.
C
C  DVAL -- (output) double precision array of values of the 
C          first derivative of the same function at the points  XVAL.
C
C  IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NVAL.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine DCHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
C                 AUGUST 1982.
C***ROUTINES CALLED  DPCHFD
C***END PROLOGUE  DPCHEV
      INTEGER  N, NVAL, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), 
     *DVAL(NVAL)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER INCFD
      LOGICAL SKIP
      DATA SKIP /.TRUE./
      DATA INCFD /1/

C
C
C***FIRST EXECUTABLE STATEMENT  DPCHEV
C
      CALL DPCHFD(N,X,F,D,INCFD,SKIP,NVAL,XVAL,FVAL,DVAL,IERR)
C
C
 5000 CONTINUE
      RETURN
      END
      SUBROUTINE DPCHEZ(N,X,F,D,SPLINE,WK,LWK,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHEZ
C***DATE WRITTEN   870821   (YYMMDD)
C***REVISION DATE  870908   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
C             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Easy to use spline or cubic Hermite interpolation.
C***DESCRIPTION
C
C          DPCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C
C     Sets derivatives for spline (two continuous derivatives) or
C     Hermite cubic (one continuous derivative) interpolation.
C     Spline interpolation is smoother, but may not "look" right if the
C     data contains both "steep" and "flat" sections.  Hermite cubics
C     can produce a "visually pleasing" and monotone interpolant to
C     monotone data. This is an easy to use driver for the routines
C     by F. N. Fritsch in reference (4) below. Various boundary
C     conditions are set to default values by DPCHEZ. Many other choices
C     are available in the subroutines PCHIC, DPCHIM and DPCHSP.
C
C     Use PCHEV to evaluate the resulting function and its derivative.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:   CALL  DPCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
C
C     INTEGER  N, IERR,  LWK
C     DOUBLE PRECISION  X(N), F(N), D(N), WK(*)
C     LOGICAL SPLINE
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(I) is value corresponding to X(I).
C
C     D -- (output) real array of derivative values at the data points.
C
C     SPLINE -- (input) logical variable to specify if the interpolant
C           is to be a spline with two continuous derivaties
C           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
C           continuous derivative (set SPLINE=.FALSE.).
C        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
C           default "not-a-knot" boundary condition, with a continuous
C           third derivative at X(2) and X(N-1). See reference (3).
C              If SPLINE=.FALSE. the interpolating Hermite cubic will be
C           monotone if the input data is monotone. Boundary conditions
C           computed from the derivative of a local quadratic unless thi
C           alters monotonicity.
C
C     WK -- (scratch) real work array, which must be declared by the cal
C           program to be at least 2*N if SPLINE is .TRUE. and not used
C           otherwise.
C
C     LWK -- (input) length of work array WK. (Error return if
C           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means tha
C                 IERR switches in the direction of monotonicity were de
C                 When SPLINE=.FALSE.,  DPCHEZ guarantees that if the inp
C                 data is monotone, the interpolant will be too. This wa
C                 is to alert you to the fact that the input data was no
C                 monotone.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
C                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
C                 PREPRINT UCRL-87559 (APRIL 1982).
C               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
C                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
C               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
C                 AUGUST 1982.
C***ROUTINES CALLED  DPCHIM,DPCHSP
C***END PROLOGUE  DPCHEZ
      INTEGER  N, LWK, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), WK(LWK)
      LOGICAL SPLINE
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER IC(2), INCFD
      DOUBLE PRECISION  VC(2)
      DATA IC(1) /0/
      DATA IC(2) /0/
      DATA INCFD /1/
C
C
C***FIRST EXECUTABLE STATEMENT  DPCHEZ
C
      IF ( SPLINE ) THEN
        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
C
C  ERROR CONDITIONS ALREADY CHECKED IN DPCHSP OR DPCHIM

      RETURN
      END
      SUBROUTINE DPCHFD(N,X,F,D,INCFD,SKIP,NE,XE,FE,DE,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHFD
C***DATE WRITTEN   811020   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHFD-S DPCHFD-D),
C             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
C             HERMITE INTERPOLATION,PIECEWISE CUBIC EVALUATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate a piecewise cubic hermite function and its first
C            derivative at an array of points.  May be used by itself
C            for Hermite interpolation, or as an evaluator for DPCHIM
C            or DPCHIC. If only function values are required, use
C            DPCHFE instead.
C***DESCRIPTION
C
C       **** Double Precision version of PCHFD ****
C
C          DPCHFD:  Piecewise Cubic Hermite Function and Derivative
C                  evaluator
C
C     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
C     gether with its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use DPCHFE, instead.
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, NE, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
C                          DE(NE)
C        LOGICAL  SKIP
C
C        CALL  DPCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
C           is the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in DPCHIM or DPCHIC).
C           SKIP will be set to .TRUE. on normal return.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real*8 array of points at which the functions are to
C           be evaluated.
C
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XE are increasing relative to X;
C              that is,   XE(J) .GE. X(I)
C              implies    XE(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XE are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C     FE -- (output) real*8 array of values of the cubic Hermite
C           function defined by  N, X, F, D  at the points  XE.
C
C     DE -- (output) real*8 array of values of the first derivative of
C           the same function at the points  XE.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NE.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine DCHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCHFDV,XERROR
C***END PROLOGUE  DPCHFD
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-03   Minor cosmetic changes for release 1.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. To produce a single precision version, simply:
C        a. Change DPCHFD to PCHFD, and DCHFDV to CHFDV, wherever they
C           occur,
C        b. Change the double precision declaration to real,
C
C     2. Most of the coding between the call to DCHFDV and the end of
C        the IR-loop could be eliminated if it were permissible to
C        assume that XE is ordered relative to X.
C
C     3. DCHFDV does not assume that X1 is less than X2.  thus, it would
C        be possible to write a version of DPCHFD that assumes a strict-
C        ly decreasing X-array by simply running the IR-loop backwards
C        (and reversing the order of appropriate tests).
C
C     4. The present code has a minor bug, which I have decided is not
C        worth the effort that would be required to fix it.
C        If XE contains points in [X(N-1),X(N)], followed by points .LT.
C        X(N-1), followed by points .GT.X(N), the extrapolation points
C        will be counted (at least) twice in the total returned in IERR.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, NE, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
     * DE(NE)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, IERC, IR, J, JFIRST, NEXT(2), NJ
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHFD
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
C
C  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
C                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
C
C     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
C
         IF (JFIRST .GT. NE)  GO TO 5000
C
C     LOCATE ALL POINTS IN INTERVAL.
C
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
C
C     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
C
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
C
   40    CONTINUE
         NJ = J - JFIRST
C
C     SKIP EVALUATION IF NO POINTS IN INTERVAL.
C
         IF (NJ .EQ. 0)  GO TO 50
C
C     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
C
C       ----------------------------------------------------------------
        CALL DCHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR)
     *              ,NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
C       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
C
         IF (NEXT(2) .EQ. 0)  GO TO 42
C        IF (NEXT(2) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
C           RIGHT OF X(IR).
C
            IF (IR .LT. N)  GO TO 41
C           IF (IR .EQ. N)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
C           ELSE
C              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
C           ENDIF
C        ENDIF
   42    CONTINUE
C
         IF (NEXT(1) .EQ. 0)  GO TO 49
C        IF (NEXT(1) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
C           LEFT OF X(IR-1).
C
            IF (IR .GT. 2)  GO TO 43
C           IF (IR .EQ. 2)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
C           ELSE
C              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
C              EVALUATION INTERVAL.
C
C              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
C              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
C                     IN DCHFDV.
               GO TO 5005
C
   45          CONTINUE
C              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
C
C              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
C              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
C
   47          CONTINUE
C              AT THIS POINT, EITHER  XE(J) .LT. X(1)
C                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
C              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
C              CYCLING.
               IR = MAX0(1, I-1)
C           ENDIF
C        ENDIF
   49    CONTINUE
C
         JFIRST = J
C
C     END OF IR-LOOP.
C
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHFD -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHFD -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -4
      CALL XERROR ('DPCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     ERROR RETURN FROM DCHFDV.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      CALL XERROR ('DPCHFD -- ERROR RETURN FROM DCHFDV -- FATAL'
     *           , 0, IERR, 2)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DPCHIA(N,X,F,D,INCFD,SKIP,A,B,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHIA
C***DATE WRITTEN   820730   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHIA-S DPCHIA-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an arbitrary interval.
C***DESCRIPTION
C
C       **** Double Precision version of PCHIA ****
C
C          DPCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [A, B].
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
C        DOUBLE PRECISION  VALUE, DPCHIA
C        LOGICAL  SKIP
C
C        VALUE = DPCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
C           is the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in DPCHIM or DPCHIC).
C           SKIP will be set to .TRUE. on return with IERR.GE.0 .
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DCHFIV,DPCHID,XERROR
C***END PROLOGUE  DPCHIA
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Corrected conversion to double precision.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHIA to PCHIA wherever it occurs,
C        b. Change DPCHID to PCHID wherever it occurs,
C        c. Change DPCHIV to PCHIV wherever it occurs,
C        d. Change the double precision declarations to real,  and
C        e. Change the constant  ZERO  to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N), A, B
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, IA, IB, IERD, IERV, IL, IR
      DOUBLE PRECISION  VALUE, XA, XB, ZERO
      DOUBLE PRECISION  DCHFIV, DPCHID
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHIA
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (A .EQ. B)  THEN
         VALUE = ZERO
      ELSE
         XA = DMIN1 (A, B)
         XB = DMAX1 (A, B)
         IF (XB .LE. X(2))  THEN
C           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
C                   --------------------------------------------
            VALUE = DCHFIV (X(1),X(2), F(1,1),F(1,2),
     *                                D(1,1),D(1,2), A, B, IERV)
C                   --------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE IF (XA .GE. X(N-1))  THEN
C           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
c
            VALUE = DCHFIV(X(N-1),X(N), F(1,N-1),F(1,N),
     *                                 D(1,N-1),D(1,N), A, B, IERV)
C                   -----------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE
C           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
C      ......LOCATE IA AND IB SUCH THAT
C               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
C             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
C             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
C
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
C             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
C             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
C
C     ......COMPUTE THE INTEGRAL.
            IERV = 0
            IF (IB .LT. IA)  THEN
C              THIS MEANS IB = IA-1 AND
C                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
C                      ------------------------------------------------
               VALUE = DCHFIV (X(IB),X(IA), F(1,IB),F(1,IA),
     *                                     D(1,IB),D(1,IA), A, B, IERV)
C                      ------------------------------------------------
               IF (IERV .LT. 0)  GO TO 5004
            ELSE
C
C              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
               IF (IB .EQ. IA)  THEN
                  VALUE = ZERO
               ELSE
C                         ---------------------------------------------
                  VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
C                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5005
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX0 (1, IA-1)
                  IR = IL + 1
C                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), XA, X(IA), IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN0 (IB+1, N)
                  IL = IR - 1
C                                 -------------------------------------
                  VALUE = VALUE + DCHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), X(IB), XB, IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
      DPCHIA = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIA -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIA -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIA -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     TROUBLE IN DCHFIV.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERROR ('DPCHIA -- TROUBLE IN DCHFIV'
     *           , 0, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     TROUBLE IN DPCHID.  (SHOULD NEVER OCCUR.)
      IERR = -5
      CALL XERROR ('DPCHIA -- TROUBLE IN DPCHID'
     *           , 0, IERR, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DPCHID(N,X,F,D,INCFD,SKIP,IA,IB,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHID
C***DATE WRITTEN   820723   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHID-S DPCHID-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an interval whose endpoints are data
C            points.
C***DESCRIPTION
C
C       **** Double Precision version of PCHID ****
C
C          DPCHID:  Piecewise Cubic Hermite Integrator, Data limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
C
C     To provide compatibility with DPCHIM and DPCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IA, IB, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
C           is the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in DPCHIM or DPCHIC).
C           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
C
C     IA,IB -- (input) indices in X-array for the limits of integration.
C           both must be in the range [1,N].  (Error return if not.)
C           No restrictions on their relative values.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IA or IB is out of range.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DPCHID
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHID to PCHID wherever it occurs,
C        b. Change the double precision declarations to real,  and
C        c. Change the constants ZERO, HALF, SIX to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, IA, IB, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, IUP, LOW
      DOUBLE PRECISION  H, HALF, SIX, SUM, VALUE, ZERO
C
C  INITIALIZE.
C
      DATA  ZERO /0.D0/,  HALF/.5D0/, SIX/6.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHID
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (IA .EQ. IB)  THEN
         VALUE = ZERO
      ELSE
         LOW = MIN0(IA, IB)
         IUP = MAX0(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +
     *                      (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
C
C  NORMAL RETURN.
C
      DPCHID = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHID -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHID -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHID -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERROR ('DPCHID -- IA OR IB OUT OF RANGE'
     *           , 0, IERR, 1)
      RETURN
      END
      SUBROUTINE DPCHIM(N,X,F,D,INCFD,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHIM
C***DATE WRITTEN   811103   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(DPCHIM-S DPCHIM-D),
C             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
C             PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Set derivatives needed to determine a monotone piecewise
C            cubic Hermite interpolant to given data.  Boundary values
C            are provided which are compatible with monotonicity.  The
C            interpolant will have an extremum at each point where mono-
C            tonicity switches direction.  (See DPCHIC if user control
C            is desired over boundary or switch conditions.)
C***DESCRIPTION
C
C       **** Double Precision version of DPCHIM ****
C
C          DPCHIM:  Piecewise Cubic Hermite Interpolation to
C                  Monotone data.
C
C     Sets derivatives needed to determine a monotone piecewise cubic
C     Hermite interpolant to the data given in X and F.
C
C     Default boundary conditions are provided which are compatible
C     with monotonicity.  (See DPCHIC if user control of boundary con-
C     ditions is desired.)
C
C     If the data are only piecewise monotonic, the interpolant will
C     have an extremum at each point where monotonicity switches direc-
C     tion.  (See DPCHIC if user control is desired in such cases.)
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by DPCHFE or DPCHFD.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C
C        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of dependent variable values to be
C           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
C           X(I).  DPCHIM is designed for monotonic data, but it will
C           work for any F-array.  It will force extrema at points where
C           monotonicity switches direction.  If some other treatment of
C           switch points is desired, DPCHIC should be used instead.
C                                     -----
C     D -- (output) real*8 array of derivative values at the data
C           points.  If the data are monotonic, these values will
C           determine a monotone cubic Hermite function.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that IERR switches in the direction
C                 of monotonicity were detected.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
C                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
C                 PREPRINT UCRL-87559 (APRIL 1982).
C***ROUTINES CALLED  DPCHST,XERROR
C***END PROLOGUE  DPCHIM
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-02-01   1. Introduced  DPCHST  to reduce possible over/under-
C                   flow problems.
C                2. Rearranged derivative formula for same reason.
C     82-06-02   1. Modified end conditions to be continuous functions
C                   of data when monotonicity switches in next interval.
C                2. Modified formulas so end conditions are less prone
C                   of over/underflow problems.
C     82-08-03   Minor cosmetic changes for release 1.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C     2. To produce a single precision version, simply:
C        a. Change DPCHIM to PCHIM wherever it occurs,
C        b. Change DPCHST to PCHST wherever it occurs,
C        c. Change all references to the Fortran intrinsics to their
C           single precision equivalents,
C        d. Change the double precision declarations to real, and
C        e. Change the constants ZERO and THREE to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER N, INCFD, IERR
      DOUBLE PRECISION  X(N), F(INCFD,N), D(INCFD,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER I, NLESS1
      DOUBLE PRECISION  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE,
     *      H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      DOUBLE PRECISION  DPCHST
      DATA  ZERO /0.D0/, THREE/3.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
C
C  NORMAL CASE  (N .GE. 3).
C
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
C
C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (DABS(D(1,1)) .GT. DABS(DMAX))  D(1,1) = DMAX
      ENDIF
C
C  LOOP THROUGH INTERIOR POINTS.
C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
C
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
C
C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
C
         D(1,I) = ZERO
         IF ( DPCHST(DEL1,DEL2) )  42, 41, 45
C
C        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
C
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( DPCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
C
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = DMAX1( DABS(DEL1), DABS(DEL2) )
         DMIN = DMIN1( DABS(DEL1), DABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
C
   50 CONTINUE
C
C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( DPCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( DPCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (DABS(D(1,N)) .GT. DABS(DMAX))  D(1,N) = DMAX
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHIM -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHIM -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DPCHQA(N,X,F,D,A,B,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHQA
C***DATE WRITTEN   870829   (YYMMDD)
C***REVISION DATE  870829   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  EASY TO USE CUBIC HERMITE OR SPLINE INTEGRATION
C             NUMERICAL INTEGRATION, QUADRATURE
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the definite integral of a piecewise cubic Hermit
C            or spline function over an arbitrary interval, easy to use.
C***DESCRIPTION
C
C          DPCHQA:  Piecewise Cubic Hermite or Spline Integrator,
C                  Arbitrary limits, Easy to Use.
C
C          From the book "Numerical Methods and Software"
C                  by  D. Kahaner, C. Moler, S. Nash
C                          Prentice Hall 1988
C
C     Evaluates the definite integral of the cubic Hermite or spline
C     function defined by  N, X, F, D  over the interval [A, B].  This
C     is an easy to use driver for the routine DPCHIA by F.N. Fritsch
C     described in reference (2) below. That routine also has other
C     capabilities.
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C           VALUE = DPCHQA (N, X, F, D, A, B, IERR)
C
C     INTEGER  N, IERR
C     DOUBLE PRECISION  X(N), F(N), D(N), A, B
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) double precision array of independent variable
C           values.  The elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) double precision array of function values.
C           F(I) is the value corresponding to X(I).
C
C     D -- (input) double precision array of derivative values.  D(I) is
C           the value corresponding to X(I).
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194
C                 AUGUST 1982.
C***ROUTINES CALLED  DPCHIA
C***END PROLOGUE  DPCHQA
      INTEGER  N, IERR
      DOUBLE PRECISION  X(N), F(N), D(N), A, B
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  INCFD
      DOUBLE PRECISION  DPCHIA
      LOGICAL SKIP
C
C  INITIALIZE.
C
      DATA  INCFD /1/
      DATA  SKIP /.TRUE./
C
C
C***FIRST EXECUTABLE STATEMENT  DPCHQA

      DPCHQA  =  DPCHIA( N, X, F, D, INCFD, SKIP, A, B, IERR )
C
C ERROR MESSAGES ARE FROM LOWER LEVEL ROUTINES
      RETURN
      END
      SUBROUTINE DPCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHSP
C***DATE WRITTEN   820503   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=DOUBLE PRECISION(PCHSP-S DPCHSP-D),
C             CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
C             SPLINE INTERPOLATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Set derivatives needed to determine the Hermite represen-
C            tation of the cubic spline interpolant to given data, with
C            specified boundary conditions.
C***DESCRIPTION
C
C       **** Double Precision version of PCHSP ****
C
C          DPCHSP:   Piecewise Cubic Hermite Spline
C
C     Computes the Hermite representation of the cubic spline inter-
C     polant to the data given in X and F satisfying the boundary
C     conditions specified by IC and VC.
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by DPCHFE or DPCHFD.
C
C     NOTE:  This is a modified version of C. de Boor'S cubic spline
C            routine CUBSPL.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  IC(2), N, NWK, IERR
C        DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
C
C        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
C
C   Parameters:
C
C     IC -- (input) integer array of length 2 specifying desired
C           boundary conditions:
C           IC(1) = IBEG, desired condition at beginning of data.
C           IC(2) = IEND, desired condition at end of data.
C
C           IBEG = 0  to set D(1) so that the third derivative is con-
C              tinuous at X(2).  This is the "not a knot" condition
C              provided by de Boor'S cubic spline routine CUBSPL.
C              < This is the default boundary condition. >
C           IBEG = 1  if first derivative at X(1) is given in VC(1).
C           IBEG = 2  if second derivative at X(1) is given in VC(1).
C           IBEG = 3  to use the 3-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.3 .)
C           IBEG = 4  to use the 4-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.4 .)
C          NOTES:
C           1. An error return is taken if IBEG is out of range.
C           2. For the "natural" boundary condition, use IBEG=2 and
C              VC(1)=0.
C
C           IEND may take on the same values as IBEG, but applied to
C           derivative at X(N).  In case IEND = 1 or 2, the value is
C           given in VC(2).
C
C          NOTES:
C           1. An error return is taken if IEND is out of range.
C           2. For the "natural" boundary condition, use IEND=2 and
C              VC(2)=0.
C
C     VC -- (input) real*8 array of length 2 specifying desired boundary
C           values, as indicated above.
C           VC(1) need be set only if IC(1) = 1 or 2 .
C           VC(2) need be set only if IC(2) = 1 or 2 .
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real*8 array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real*8 array of dependent variable values to be
C           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
C           X(I).
C
C     D -- (output) real*8 array of derivative values at the data
C           points.  These values will determine the cubic spline
C           interpolant with the requested boundary conditions.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     WK -- (scratch) real*8 array of working storage.
C
C     NWK -- (input) length of work array.
C           (Error return if NWK.LT.2*N .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
C              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
C              IERR = -6  if both of the above are true.
C              IERR = -7  if NWK is too small.
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C             (The D-array has not been changed in any of these cases.)
C              IERR = -8  in case of trouble solving the linear system
C                         for the interior derivative values.
C             (The D-array may have been changed in this case.)
C             (             Do **NOT** use it!                )
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
C                 VERLAG (NEW YORK, 1978), PP. 53-59.
C***ROUTINES CALLED  DPCHDF,XERROR
C***END PROLOGUE  DPCHSP
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Corrected XERROR calls for d.p. name(s).
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHSP to PCHSP wherever it occurs,
C        b. Change the double precision declarations to real, and
C        c. Change the constants ZERO, HALF, ... to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER IC(2), N, INCFD, NWK, IERR
      DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER IBEG, IEND, INDEX, J, NM1
      DOUBLE PRECISION  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4),
     *  ZERO
      DOUBLE PRECISION  DPCHDF
C
      DATA  ZERO /0.D0/, HALF/.5D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DPCHSP
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
      IF ( NWK .LT. 2*N )  GO TO 5007
C
C  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
C  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
C
C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
C
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
C
C  SET UP FOR BOUNDARY CONDITIONS.
C
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
C        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
C           INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
C                 --------------------------------
         D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
C
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
C        PICK UP LAST IEND POINTS.
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
C           INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
C                 --------------------------------
         D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
C
C --------------------( BEGIN CODING FROM CUBSPL )--------------------
C
C  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
C  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
C  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
C     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
C
C  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
C             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
C
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
C           NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
C           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)
     *                        + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
C        SLOPE PRESCRIBED AT LEFT END.
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
C
C  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
C  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
C  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
C
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1)
     *                  + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
C
C  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
C           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
C
C     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
C     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
C     AT THIS POINT.
      IF (IEND .EQ. 1)  GO TO 30
C
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
C           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
C           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
C           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
C           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
C           KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
C           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     *                  + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
C
C  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
C
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
C
C  CARRY OUT BACK SUBSTITUTION
C
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
C --------------------(  END  CODING FROM CUBSPL )--------------------
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('DPCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 0, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('DPCHSP -- INCREMENT LESS THAN ONE'
     *           , 0, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('DPCHSP -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 0, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERROR ('DPCHSP -- IC OUT OF RANGE'
     *           , 0, IERR, 1)
      RETURN
C
 5007 CONTINUE
C     NWK TOO SMALL RETURN.
      IERR = -7
      CALL XERROR ('DPCHSP -- WORK ARRAY TOO SMALL'
     *           , 0, IERR, 1)
      RETURN
C
 5008 CONTINUE
C     SINGULAR SYSTEM.
C   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
C   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      CALL XERROR ('DPCHSP -- SINGULAR LINEAR SYSTEM'
     *           , 0, IERR, 1)
      RETURN
C
 5009 CONTINUE
C     ERROR RETURN FROM DPCHDF.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERROR ('DPCHSP -- ERROR RETURN FROM DPCHDF'
     *           , 0, IERR, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DPCHST(ARG1,ARG2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DPCHST
C***REFER TO  DPCHCE,DPCHCI,DPCHCS,DPCHIM
C***ROUTINES CALLED  (NONE)
C***DESCRIPTION
C
C         DPCHST:  DPCHIP Sign-Testing Routine.
C
C
C     Returns:
C        -1. if ARG1 and ARG2 are of opposite sign.
C         0. if either argument is zero.
C        +1. if ARG1 and ARG2 are of the same sign.
C
C     The object is to do this without multiplying ARG1*ARG2, to avoid
C     possible over/underflow problems.
C
C  Fortran intrinsics used:  SIGN.
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DPCHST to PCHST wherever it occurs,
C        b. Change all references to the Fortran intrinsics to their
C           single presision equivalents,
C        c. Change the double precision declarations to real, and
C        d. Change the constants  ZERO  and  ONE  to single precision.
C***END PROLOGUE  DPCHST
C
C  DECLARE ARGUMENTS.
C
      DOUBLE PRECISION  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      DOUBLE PRECISION  ONE, ZERO
      DATA  ZERO /0.D0/,  ONE/1.D0/
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  DPCHST
      DPCHST = DSIGN(ONE,ARG1) * DSIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  DPCHST = ZERO

      RETURN
      END
      SUBROUTINE DQ1DA(f,A,B,EPS,R,E,KF,IFLAG)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQ1DA
C***DATE WRITTEN  870525
C***REVISION DATE  870525
C***CATEGORY NO.  H2A1A1
C***KEYWORDS   ADAPTIVE  QUADRATURE, AUTOMATIC  QUADRATURE
C***AUTHOR  KAHANER, DAVID K., SCIENTIFIC COMPUTING DIVISION, NBS.
C***PURPOSE  APPROXIMATES ONE DIMENSIONAL INTEGRALS OF USER DEFINED
C            FUNCTIONS, EASY TO USE. DOUBLE PRECISION VERSION OF Q1DA.
C
C***DESCRIPTION
C       DQ1DA IS A SUBROUTINE FOR THE AUTOMATIC EVALUATION
C            OF THE DEFINITE INTEGRAL OF A USER DEFINED FUNCTION
C            OF ONE VARIABLE. DOUBLE PRECISION VERSION OF Q1DA
C
C       A R G U M E N T S   I N   T H E   C A L L   S E Q U E N C E
C
C       F     (INPUT), the name of the function that evaluates
C             the integrand.
C       A
C       B     (INPUT) THE ENDPOINTS OF THE INTEGRATION INTERVAL
C       EPS   (INPUT) THE ACCURACY TO WHICH YOU WANT THE INTEGRAL
C                COMPUTED.  IF YOU WANT 2 DIGITS OF ACCURACY SET
C                EPS=.01, FOR 3 DIGITS SET EPS=.001, ETC.
C                EPS MUST BE POSITIVE.
C       R     (OUTPUT) DQ1DA'S BEST ESTIMATE OF YOUR INTEGRAL
C       E     (OUTPUT) AN ESTIMATE OF ABS(INTEGRAL-R)
C       KF    (OUTPUT) THE COST OF THE INTEGRATION, MEASURED IN
C                   NUMBER OF EVALUATIONS OF YOUR INTEGRAND.
C                   KF WILL ALWAYS BE AT LEAST 30.
C       IFLAG (OUTPUT) TERMINATION FLAG...POSSIBLE VALUES ARE
C               0   NORMAL COMPLETION, E SATISFIES
C                        E<EPS  AND  E<EPS*ABS(R)
C               1   NORMAL COMPLETION, E SATISFIES
C                        E<EPS, BUT E>EPS*ABS(R)
C               2   NORMAL COMPLETION, E SATISFIES
C                        E<EPS*ABS(R), BUT E>EPS
C               3   NORMAL COMPLETION BUT EPS WAS TOO SMALL TO
C                     SATISFY ABSOLUTE OR RELATIVE ERROR REQUEST.
C
C               4   ABORTED CALCULATION BECAUSE OF SERIOUS ROUNDING
C                     ERROR.  PROBABLY E AND R ARE CONSISTENT.
C               5   ABORTED CALCULATION BECAUSE OF INSUFFICIENT STORAGE.
C                     R AND E ARE CONSISTENT.
C               6   ABORTED CALCULATION BECAUSE OF SERIOUS DIFFICULTIES
C                     MEETING YOUR ERROR REQUEST.
C               7   ABORTED CALCULATION BECAUSE EPS WAS SET <= 0.0
C
C            NOTE...IF IFLAG=3, 4, 5 OR 6 CONSIDER USING DQ1DAX INSTEAD.
C
C       NOTE... A,B,EPS, R AND E MUST BE DECLARED DOUBLE PRECISION IN
C                YOUR CALLING PROGRAM.
C    W H E R E   I S   Y O U R   I N T E G R A N D ?
C
C        YOU MUST WRITE A FORTRAN FUNCTION, CALLED F, TO EVALUATE
C        THE INTEGRAND.  USUALLY THIS LOOKS LIKE...
C                 DOUBLE PRECISION FUNCTION F(X)
C                    DOUBLE PRECISION X
C                    F=(EVALUATE THE INTEGRAND AT THE POINT X)
C                    RETURN
C                 END
C
C
C    T Y P I C A L   P R O B L E M   S E T U P
C
C          DOUBLE PRECISION A,B,EPS,R,E,F
C          EXTERNAL F
C
C          A=0.0
C          B=1.0          (SET INTERVAL ENDPOINTS TO [0,1])
C          EPS=0.001       (SET ACCURACY REQUEST FOR 3 DIGITS)
C          CALL DQ1DA(F,A,B,EPS,R,E,KF,IFLAG)
C          END
C          DOUBLE PRECISION FUNCTION F(X)
C            DOUBLE PRECISION X
C            F=SIN(2.*X)-SQRT(X)     (FOR EXAMPLE)
C            RETURN
C          END
C      FOR THIS SAMPLE PROBLEM, THE OUTPUT IS
C  0.0    1.0     .001    .041406750    .69077D-07    30    0
C
C    R E M A R K   I.
C
C           A SMALL AMOUT OF RANDOMIZATION IS BUILT INTO THIS PROGRAM.
C           CALLING DQ1DA A FEW TIMES IN SUCCESSION WILL GIVE DIFFERENT
C           BUT HOPEFULLY CONSISTENT RESULTS.
C
C   R E M A R K   II.
C
C           THIS ROUTINE IS DESIGNED FOR INTEGRATION OVER A FINITE
C           INTERVAL.  THUS THE INPUT ARGUMENTS A AND B MUST BE
C           VALID REAL NUMBERS ON YOUR COMPUTER.  IF YOU WANT TO DO
C           AN INTEGRAL OVER AN INFINITE INTERVAL SET A OR B OR BOTH
C           LARGE ENOUGH SO THAT THE INTERVAL [A,B] CONTAINS MOST OF
C           THE INTEGRAND.  CARE IS NECESSARY, HOWEVER.  FOR EXAMPLE,
C           TO INTEGRATE EXP(-X*X) ON THE ENTIRE REAL LINE ONE COULD
C           TAKE A=-20., B=20. OR SIMILAR VALUES TO GET GOOD RESULTS.
C           IF YOU TOOK A=-1.E10 AND B=+1.E10 TWO BAD THINGS WOULD
C           OCCUR. FIRST, YOU WILL CERTAINLY GET AN ERROR MESSAGE FROM
C           THE EXP ROUTINE, AS ITS ARGUMENT IS TOO SMALL.  OTHER
C           THINGS COULD HAPPEN TOO, FOR EXAMPLE AN UNDERFLOW.
C           SECOND, EVEN IF THE ARITHMETIC WORKED PROPERLY DQ1DA WILL
C           SURELY GIVE AN INCORRECT ANSWER, BECAUSE ITS FIRST TRY
C           AT SAMPLING THE INTEGRAND IS BASED ON YOUR SCALING AND
C           IT IS VERY UNLIKELY TO SELECT EVALUATION POINTS IN THE
C           INFINITESMALLY SMALL INTERVAL [-20,20] WHERE ALL THE
C           INTEGRAND IS CONCENTRATED, WHEN A, B ARE SO LARGE.
C
C    M O R E   F L E X I B I L I T Y
C
C           DQ1DA IS AN EASY TO USE DRIVER FOR ANOTHER PROGRAM, DQ1DAX.
C           DQ1DAX PROVIDES SEVERAL OPTIONS WHICH ARE NOT AVAILABLE
C                WITH DQ1DA.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED (DQ1DAX)
C    E N D   O F   D O C U M E N T A T I O N
C***END PROLOGUE  DQ1DA
C
      DOUBLE PRECISION A,B,E,EPS
      external f
      double precision f
      double precision FMAX,FMIN,R,W(50,6)
      LOGICAL RST
C
C***FIRST EXECUTUTABLE STATEMENT DQ1DA

      NINT=1
      RST = .FALSE.
      NMAX=50
      CALL DQ1DAX(F,A,B,EPS,R,E,NINT,RST,W,NMAX,FMIN,FMAX,KF,IFLAG)

      RETURN
      END
      SUBROUTINE DQ1DAX(
     *      F,A,B,EPS,R,E,NINT,RST,W,NMAX,FMIN,FMAX,KF,IFLAG)

c*********************************************************************72
c
C***BEGIN PROLOGUE DQ1DAX
C***DATE WRITTEN 870525
C***REVISION DATE 870525
C***CATEGORY NO. H2A1A1
C***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC QUADRATURE
C***AUTHOR  KAHANER, DAVID K., SCIENTIFIC COMPUING DIVISION, NBS.
C***PURPOSE  APPROXIMATES ONE DIMENSIONAL INTEGRALS OF USER DEFINED
C            FUNCTIONS, FLEXIBLE USAGE. DOUBLE PRECISION VERSION OF
C            Q1DAX.
C
C***DESCRIPTION
C
C       DQ1DAX IS A FLEXIBLE SUBROUTINE FOR THE AUTOMATIC EVALUATION
C             OF DEFINITE INTEGRALS OF A USER DEFINED FUNCTION
C             OF ONE VARIABLE. DOUBLE PRECISION VERSION OF Q1DAX.
C
C             FOR AN EASIER TO USE ROUTINE SEE DQ1DA.
C
C             CAPABILITIES OF DQ1DAX (IN ADDITION TO THOSE OF DQ1DA)
C             INCLUDE:
C                ABILITY TO RESTART A CALCULATION TO GREATER
C                  ACCURACY WITHOUT PENALTY...
C                ABILITY TO SPECIFY AN INITIAL PARTITION OF
C                  THE INTEGRATION INTERVAL...
C                ABILITY TO INCREASE THE WORK SPACE TO HANDLE
C                  MORE DIFFICULT PROBLEMS...
C                OUTPUT OF LARGEST/SMALLEST INTEGRAND VALUE FOR
C                  APPLICATIONS SUCH AS SCALING GRAPHS...
C
C       A R G U M E N T S   I N   T H E   C A L L   S E Q U E N C E
C
C       F      (INPUT) THE NAME OF YOUR INTEGRAND FUNCTION.
C                      THIS NAME MUST APPEAR IN AN EXTERNAL STATEMENT
C                      IN ANY PROGRAM WHICH CALLS DQ1DAX.
C                      YOU MUST WRITE F IN THE FORM
C                           DOUBLE PRECISION FUNCTION F(X)
C                              DOUBLE PRECISION X
C                              F=(EVALUATE INTEGRAND AT THE POINT X)
C                              RETURN
C                           END
C       A
C       B       (INPUT) ENDPOINTS OF INTEGRATION INTERVAL
C       EPS     (INPUT)  ACCURACY TO WHICH THE INTEGRAL IS TO BE CALCULATED.
C                          DQ1DAX WILL TRY TO ACHIEVE RELATIVE ACCURACY,
C                          E.G. SET EPS=.01 FOR 2 DIGITS, .001 FOR 3, ETC.
C       R       (OUTPUT) THE ESTIMATE OF THE INTEGRAL
C       E       (OUTPUT) THE ESTIMATE OF THE ABSOLUTE ERROR IN R.
C       NINT    (INPUT
C                OUTPUT)
C                        AS AN INPUT QUANTITY, NINT MUST BE SET TO
C                          THE NUMBER OF SUBINTERVALS IN THE INITIAL
C                          PARTITION OF [A,B].  FOR MOST PROBLEMS
C                          THIS IS JUST 1, THE INTERVAL [A,B] ITSELF.
C                          NINT MUST BE LESS THAN NMAX, SEE BELOW.
C                          NINT IS USEFUL IF YOU WOULD LIKE TO HELP
C                          DQ1DAX LOCATE A DIFFICULT SPOT ON [A,B].
C                          IN THIS REGARD NINT IS USED ALONG
C                          WITH THE ARRAY W (SEE BELOW).  IF YOU SET
C                          NINT=1 IT IS NOT NECESSARY TO BE CONCERNED
C                          WITH W, EXCEPT THAT IT MUST BE DIMENSIONED...
C                          AS AN EXAMPLE OF MORE GENERAL APPLICATIONS,
C                          IF [A,B]=[0,1] BUT THE INTEGRAND JUMPS AT 0.3,
C                          IT WOULD BE WISE TO SET NINT=2 AND THEN SET
C                                  W(1,1)=0.0  (LEFT ENDPOINT)
C                                  W(2,1)=0.3  (SINGULAR POINT)
C                                  W(3,1)=1.0  (RIGHT ENDPOINT)
C                          IF YOU SET NINT GREATER THAN 1, BE SURE TO
C                          CHECK THAT YOU HAVE ALSO SET
C                                W(1,1)=A  AND  W(NINT+1,1)=B
C                        AS AN OUTPUT QUANTITY, NINT GIVES THE
C                          NUMBER OF SUBINTERVALS IN THE FINAL
C                          PARTITION OF [A,B].
C       RST     (INPUT) A LOGICAL VARIABLE (E.G. TRUE OR FALSE)
C                       SET RST=.FALSE. FOR INITIAL CALL TO DQ1DAX
C                       SET RST=.TRUE. FOR A SUBSEQUENT CALL,
C                           E.G. ONE FOR WHICH MORE ACCURACY IS
C                           DESIRED (SMALLER EPS).  A RESTART ONLY
C                           MAKES SENSE IF THE PRECEDING CALL RETURNED
C                           WITH A VALUE OF IFLAG (SEE BELOW) LESS THAN 3.
C                           ON A RESTART YOU MAY NOT CHANGE THE VALUES OF ANY
C                           OTHER ARGUMENTS IN THE CALL SEQUENCE, EXCEPT EPS.
C       W(NMAX,6) (INPUT &
C                  SCRATCH)
C                       W IS AN ARRAY USED BY DQ1DAX.
C                        YOU  M U S T  INCLUDE A DIMENSION STATEMENT IN
C                        YOUR CALLING PROGRAM TO ALLOCATE THIS STORAGE.
C                        THIS SHOULD BE OF THE FORM
C                                   DIMENSION W(NMAX,6)
C                        WHERE NMAX IS AN INTEGER. AN ADEQUATE VALUE OF
C                        NMAX IS 50.  IF YOU SET NINT>1 YOU MUST ALSO
C                        INITIALIZE W, SEE NINT ABOVE.
C       NMAX    (INPUT) AN INTEGER EQUAL TO THE FIRST SUBSCRIPT IN THE
C                        DIMENSION STATEMENT FOR THE ARRAY W.  THIS IS
C                        ALSO EQUAL TO THE MAXIMUM NUMBER OF SUBINTERVALS
C                        PERMITTED IN THE INTERNAL PARTITION OF [A,B].
C                        A VALUE OF 50 IS AMPLE FOR MOST PROBLEMS.
C       FMIN
C       FMAX    (OUTPUT) THE SMALLEST AND LARGEST VALUES OF THE INTEGRAND
C                          WHICH OCCURRED DURING THE CALCULATION.  THE
C                          ACTUAL INTEGRAND RANGE ON [A,B] MAY, OF COURSE,
C                          BE GREATER BUT PROBABLY NOT BY MORE THAN 10%.
C       KF      (OUTPUT) THE ACTUAL NUMBER OF INTEGRAND EVALUATIONS USED
C                          BY DQ1DAX TO APPROXIMATE THIS INTEGRAL.  KF
C                          WILL ALWAYS BE AT LEAST 30.
C       IFLAG   (OUTPUT) TERMINATION FLAG...POSSIBLE VALUES ARE
C                  0   NORMAL COMPLETION, E SATISFIES
C                           E<EPS  AND  E<EPS*ABS(R)
C                  1   NORMAL COMPLETION, E SATISFIES
C                           E<EPS, BUT E>EPS*ABS(R)
C                  2   NORMAL COMPLETION, E SATISFIES
C                           E<EPS*ABS(R), BUT E>EPS
C                  3   NORMAL COMPLETION BUT EPS WAS TOO SMALL TO
C                        SATISFY ABSOLUTE OR RELATIVE ERROR REQUEST.
C                  4   ABORTED CALCULATION BECAUSE OF SERIOUS ROUNDING
C                        ERROR.  PROBABLY E AND R ARE CONSISTENT.
C                  5   ABORTED CALCULATION BECAUSE OF INSUFFICIENT STORAGE.
C                        R AND E ARE CONSISTENT.  PERHAPS INCREASING NMAX
C                        WILL PRODUCE BETTER RESULTS.
C                  6   ABORTED CALCULATION BECAUSE OF SERIOUS DIFFICULTIES
C                        MEETING YOUR ERROR REQUEST.
C                  7   ABORTED CALCULATION BECAUSE EITHER EPS, NINT OR NMAX
C                        HAS BEEN SET TO AN ILLEGAL VALUE.
C                  8   ABORTED CALCULATION BECAUSE YOU SET NINT>1 BUT FORGOT
C                        TO SET W(1,1)=A  AND  W(NINT+1,1)=B
C       NOTE... A,B,EPS,R,E,W,FMIN AND FMAX MUST BE DECLARED DOUBLE
C                  PRECISION IN YOUR CALLING PROGRAM.
C
C     T Y P I C A L   P R O B L E M   S E T   U P
C
C      DOUBLE PRECISION A,B,EPS,W(50,6),R,E,FMIN,FMAX
C      LOGICAL RST
C      EXTERNAL F
C      A=0.0
C      B=1.0
C      W(1,1)=A
C      W(2,1)=.3      [SET INTERNAL PARTITION POINT AT .3]
C      W(3,1)=B
C      NINT=2         [INITIAL PARTITION HAS 2 INTERVALS]
C      RST=.FALSE.
C      EPS=.001
C      NMAX=50
C
C    1 CALL DQ1DAX(F,A,B,EPS,R,E,NINT,RST,W,NMAX,FMIN,FMAX,KF,IFLAG)
C
C      IF(EPS.EQ. .0001 .OR. IFLAG.GE.3)STOP
C      RST=.TRUE
C      EPS=.0001      [ASK FOR ANOTHER DIGIT]
C      GO TO 1
C      END
C      DOUBLE PRECISION FUNCTION F(X)
C      DOUBLE PRECISION X
C      IF(X.LT. .3)
C     1  THEN
C          F=X**(0.2)*ALOG(X)
C        ELSE
C          F=SIN(X)
C      ENDIF
C      RETURN
C      END
C
C
C            R E M A R K
C
C               WHEN YOU USE DQ1ADX WITH NINT=1, WE HAVE BUILT A SMALL
C            AMOUNT OF RANDOMIZATION INTO IT.  REPEATED CALLS DURING
C            THE SAME RUN WILL PRODUCE DIFFERENT, BUT HOPEFULLY
C            CONSISTENT, RESULTS.
C
C    E N D   O F   D O C U M E N T A T I O N
C***REFERENCES  (NONE)
C***ROUTINES CALLED (D1MACH,UNI,DGL15T,IDAMAX)
C***END PROLOGUE DQ1DAX
C
      DOUBLE PRECISION A,B,E,EB,EPMACH,EPS,F,FMAX,FMAXL,FMAXR,FMIN,
     *  FMINL,FMINR,FMN,FMX,R,D1MACH,RAB,RABS,RAV,T,TE,TE1,TE2,TR,
     *  TR1,TR2,UFLOW,W(NMAX,6),XM
      INTEGER C
      EXTERNAL F
      LOGICAL RST
C
C***FIRST EXECUTABLE STATEMENT  DQ1DAX
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
      MXTRY=NMAX/2
C          In case there is no more room, we can toss out easy intervals,
C             at most MXTRY times.
      IF(A.EQ.B) THEN
          R=0.
          E=0.
          NINT=0
          IFLAG=0
          KF=1
          FMIN=F(A)
          FMAX=FMIN
          GO TO 20
      ENDIF
      IF(RST) THEN
         IF(IFLAG.LT.3) THEN
           EB=MAX(100.*UFLOW,MAX(EPS,50.*EPMACH)*ABS(R))
           DO 19 I=1,NINT
               IF(ABS(W(I,3)).GT.(EB*(W(I,2)-W(I,1))/(B-A)))THEN
                                 W(I,3)=ABS(W(I,3))
               ELSE
                                 W(I,3)=-ABS(W(I,3))
               ENDIF
   19      CONTINUE
           GOTO 15
         ELSE
           GOTO 20
         ENDIF
      ENDIF
      KF=0
      IF(EPS .LE. 0. .OR. NINT .LE. 0 .OR. NINT .GE. NMAX) THEN
          IFLAG=7
          GO TO 20
      ENDIF
      IF(NINT.EQ.1)
     1  THEN
          W(1,1)=A
          W(2,2)=B
          W(1,5)=A
          W(1,6)=B
          W(2,5)=A
          W(2,6)=B
C          SELECT FIRST SUBDIVISION RANDOMLY
          W(1,2)=A+(B-A)/2.*(2*DUNI()+7.)/8.
          W(2,1)=W(1,2)
          NINT=2
        ELSE
          IF(W(1,1).NE.A .OR. W(NINT+1,1).NE.B) THEN
               IFLAG=8
               GO TO 20
          ENDIF
          W(1,5)=A
          DO 89 I=1,NINT
             W(I,2)=W(I+1,1)
             W(I,5)=W(I,1)
             W(I,6)=W(I,2)
   89     CONTINUE
      ENDIF
C
      IFLAG = 0
      IROFF=0
      RABS=0.0
      DO 3 I=1,NINT
          CALL DGL15T(F,W(I,1),W(I,2),W(I,5),W(I,6),
     1          W(I,4),W(I,3),RAB,RAV,FMN,FMX)
          KF=KF+15
          IF(I.EQ.1)
     1       THEN
               R=W(I,4)
               E=W(I,3)
               RABS=RABS+RAB
               FMIN=FMN
               FMAX=FMX
             ELSE
               R=R+W(I,4)
               E=E+W(I,3)
               RABS=RABS+RAB
               FMAX=MAX(FMAX,FMX)
               FMIN=MIN(FMIN,FMN)
          ENDIF
    3 CONTINUE
      DO 10 I=NINT+1,NMAX
          W(I,3) = 0.
   10 CONTINUE
   15 CONTINUE
C
C   MAIN SUBPROGRAM LOOP
C
      IF(100.*EPMACH*RABS.GE.ABS(R) .AND. E.LT.EPS)GO TO 20
      EB=MAX(100.*UFLOW,MAX(EPS,50.*EPMACH)*ABS(R))
      IF(E.LE.EB) GO TO 20
      IF (NINT.LT.NMAX)
     1 THEN
        NINT = NINT+1
        C = NINT
       ELSE
        C=0
   16   IF(C.EQ.NMAX .OR. MXTRY.LE.0) THEN
            IFLAG=5
            GO TO 20
        ENDIF
        C=C+1
        IF(W(C,3).GT.0.0) GO TO 16
C            Found an interval to throw out
        MXTRY=MXTRY-1
      END IF
      LOC=IDAMAX(NINT,W(1,3),1)
      XM = W(LOC,1)+(W(LOC,2)-W(LOC,1))/2.
      IF ((MAX(ABS(W(LOC,1)),ABS(W(LOC,2)))).GT.
     1   ((1.+100.*EPMACH)*(ABS(XM)+0.1D+04*UFLOW)))
     2    THEN
            CALL DGL15T(F,W(LOC,1),XM,W(LOC,5),W(LOC,6),
     1                    TR1,TE1,RAB,RAV,FMINL,FMAXL)
            KF=KF+15
            IF (TE1.LT.(EB*(XM-W(LOC,1))/(B-A))) TE1=-TE1
            CALL DGL15T(F,XM,W(LOC,2),W(LOC,5),W(LOC,6),
     1                    TR2,TE2,RAB,RAV,FMINR,FMAXR)
            KF=KF+15
            FMIN=MIN(FMIN,FMINL,FMINR)
            FMAX=MAX(FMAX,FMAXL,FMAXR)
            IF (TE2.LT.(EB*(W(LOC,2)-XM)/(B-A))) TE2=-TE2
            TE = ABS(W(LOC,3))
            TR = W(LOC,4)
            W(C,3) = TE2
            W(C,4) = TR2
            W(C,1) = XM
            W(C,2) = W(LOC,2)
            W(C,5) = W(LOC,5)
            W(C,6) = W(LOC,6)
            W(LOC,3) = TE1
            W(LOC,4) = TR1
            W(LOC,2) = XM
            E = E-TE+(ABS(TE1)+ABS(TE2))
            R = R-TR+(TR1+TR2)
            IF(ABS(ABS(TE1)+ABS(TE2)-TE).LT.0.001*TE) THEN
                IROFF=IROFF+1
                IF(IROFF.GE.10) THEN
                     IFLAG=4
                     GO TO 20
                ENDIF
            ENDIF
          ELSE
            IF (EB.GT.W(LOC,3))
     1         THEN
                   W(LOC,3) = 0.
               ELSE
                   IFLAG=6
                   GO TO 20
            END IF
      END IF
      GO TO 15
C
C        ALL EXITS FROM HERE
C
   20 CONTINUE
      IF(IFLAG.GE.4)RETURN
      IFLAG=3
      T=EPS*ABS(R)
      IF(E.GT.EPS .AND. E.GT.T)RETURN
      IFLAG=2
      IF(E.GT.EPS .AND. E.LT.T)RETURN
      IFLAG=1
      IF(E.LT.EPS .AND. E.GT.T)RETURN
      IFLAG=0
      RETURN
      END
      SUBROUTINE DQAGI(F,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     1   IER,LIMIT,LENW,LAST,IWORK,WORK)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQAGI
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A3A1,H2A4A1
C***KEYWORDS  AUTOMATIC INTEGRATOR,EXTRAPOLATION,GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE,INFINITE INTERVALS,TRANSFORMATION
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            INTEGRAL   I = Integral of F over (BOUND,+INFINITY)
C            OR I = Integral of F over (-INFINITY,BOUND)
C            OR I = Integral of F over (-INFINITY,+INFINITY)
C            Hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***DESCRIPTION
C
C        Integration over infinite intervals
C        Standard fortran subroutine
C
C        PARAMETERS
C         ON ENTRY
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            BOUND  - Double precision
C                     Finite bound of integration range
C                     (has no meaning if interval is doubly-infinite)
C
C            INF    - Integer
C                     indicating the kind of integration range involved
C                     INF = 1 corresponds to  (BOUND,+INFINITY),
C                     INF = -1            to  (-INFINITY,BOUND),
C                     INF = 2             to (-INFINITY,+INFINITY).
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                   - IER.GT.0 abnormal termination of the routine. The
C                             estimates for result and error are less
C                             reliable. It is assumed that the requested
C                             accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However, if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties. If
C                             the position of a local difficulty can be
C                             determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table.
C                             It is assumed that the requested tolerance
C                             cannot be achieved, and that the returned
C                             RESULT is the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                              or LIMIT.LT.1 or LENIW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set to
C                             zero. Exept when LIMIT or LENIW is
C                             invalid, IWORK(1), WORK(LIMIT*2+1) and
C                             WORK(LIMIT*3+1) are set to ZERO, WORK(1)
C                             is set to A and WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    Dimensioning parameter for IWORK
C                    LIMIT determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    If LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for WORK
C                    LENW must be at least LIMIT*4.
C                    If LENW.LT.LIMIT*4, the routine will end
C                    with IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least LIMIT, the first
C                    K elements of which contain pointers
C                    to the error estimates over the subintervals,
C                    such that WORK(LIMIT*3+IWORK(1)),... ,
C                    WORK(LIMIT*3+IWORK(K)) form a decreasing
C                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
C                    K = LIMIT+1-LAST otherwise
C
C            WORK  - Double precision
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left
C                     end points of the subintervals in the
C                     partition of (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) Contain
C                     the right end points,
C                    WORK(LIMIT*2+1), ...,WORK(LIMIT*2+LAST) contain the
C                     integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3)
C                     contain the error estimates.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DQAGIE,XERROR
C***END PROLOGUE  DQAGI
C
      DOUBLE PRECISION ABSERR,BOUND,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,INF,IWORK,LAST,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(LIMIT),WORK(LENW)
C
      EXTERNAL F
C
C         CHECK VALIDITY OF LIMIT AND LENW.
C
C***FIRST EXECUTABLE STATEMENT  DQAGI
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF(LIMIT.LT.1.OR.LENW.LT.LIMIT*4) GO TO 10
C
C         PREPARE CALL FOR DQAGIE.
C
      L1 = LIMIT+1
      L2 = LIMIT+L1
      L3 = LIMIT+L2
C
      CALL DQAGIE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1  NEVAL,IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C         CALL ERROR HANDLER IF NECESSARY.
C
       LVL = 0
10    IF(IER.EQ.6) LVL = 1
      IF(IER.NE.0) CALL XERROR( 'ABNORMAL RETURN FROM DQAGI',26,IER,LVL)
      RETURN
      END
      SUBROUTINE DQAGIE(F,BOUND,INF,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,
     1   NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQAGIE
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A3A1,H2A4A1
C***KEYWORDS  AUTOMATIC INTEGRATOR,EXTRAPOLATION,GENERAL-PURPOSE,
C             GLOBALLY ADAPTIVE,INFINITE INTERVALS,TRANSFORMATION
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The routine calculates an approximation result to a given
C            integral   I = Integral of F over (BOUND,+INFINITY)
C            or I = Integral of F over (-INFINITY,BOUND)
C            or I = Integral of F over (-INFINITY,+INFINITY),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I))
C***DESCRIPTION
C
C Integration over infinite intervals
C Standard fortran subroutine
C
C            F      - Double precision
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            BOUND  - Double precision
C                     Finite bound of integration range
C                     (has no meaning if interval is doubly-infinite)
C
C            INF    - Double precision
C                     Indicating the kind of integration range involved
C                     INF = 1 corresponds to  (BOUND,+INFINITY),
C                     INF = -1            to  (-INFINITY,BOUND),
C                     INF = 2             to (-INFINITY,+INFINITY).
C
C            EPSABS - Double precision
C                     Absolute accuracy requested
C            EPSREL - Double precision
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.1
C
C         ON RETURN
C            RESULT - Double precision
C                     Approximation to the integral
C
C            ABSERR - Double precision
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                   - IER.GT.0 Abnormal termination of the routine. The
C                             estimates for result and error are less
C                             reliable. It is assumed that the requested
C                             accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). However,if
C                             this yields no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (e.g. SINGULARITY,
C                             DISCONTINUITY within the interval) one
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             integrator on the subranges. If possible,
C                             an appropriate special-purpose integrator
C                             should be used, which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                             The error may be under-estimated.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 4 The algorithm does not converge.
C                             Roundoff error is detected in the
C                             extrapolation table.
C                             It is assumed that the requested tolerance
C                             cannot be achieved, and that the returned
C                             result is the best which can be obtained.
C                         = 5 The integral is probably divergent, or
C                             slowly convergent. It must be noted that
C                             divergence can occur with any other value
C                             of IER.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
C                             ELIST(1) and IORD(1) are set to zero.
C                             ALIST(1) and BLIST(1) are set to 0
C                             and 1 respectively.
C
C            ALIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the left
C                     end points of the subintervals in the partition
C                     of the transformed integration range (0,1).
C
C            BLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the right
C                     end points of the subintervals in the partition
C                     of the transformed integration range (0,1).
C
C            RLIST  - Double precision
C                     Vector of dimension at least LIMIT, the first
C                      LAST  elements of which are the integral
C                     approximations on the subintervals
C
C            ELIST  - Double precision
C                     Vector of dimension at least LIMIT,  the first
C                     LAST elements of which are the moduli of the
C                     absolute error estimates on the subintervals
C
C            IORD   - Integer
C                     Vector of dimension LIMIT, the first K
C                     elements of which are pointers to the
C                     error estimates over the subintervals,
C                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
C                     form a decreasing sequence, with K = LAST
C                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
C                     otherwise
C
C            LAST   - Integer
C                     Number of subintervals actually produced
C                     in the subdivision process
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DEA,DQK15I,DQPSRT,D1MACH
C***END PROLOGUE  DQAGIE
      DOUBLE PRECISION ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,
     1  A2,BLIST,BOUN,BOUND,B1,B2,CORREC,DABS,DEFABS,DEFAB1,DEFAB2,
     2  DMAX1,DRES,D1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,
     3  ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,ERTEST,F,RESABS,
     4  RESEPS,RESULT,RLIST,RLIST2,SMALL,UFLOW
      INTEGER ID,IER,IERR,IERRO,INF,IORD,IROFF1,IROFF2,IROFF3,JUPBND,K,
     1  KSGN,KTMIN,LAST,LIMEXP,LIMIT,MAXERR,NEVAL,NRMAX
      LOGICAL EXTRAP,LERR,NEWFLG,NOEXT
C
      PARAMETER (LIMEXP = 50)
C
C           LIMEXP IS THE SIZE OF THE EPSILON TABLE THAT CAN BE
C           GENERATED IN EA
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     1  RLIST(LIMIT),RLIST2(LIMEXP+7)
C
      EXTERNAL F
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           RLIST2    - ARRAY OF DIMENSION AT LEAST (LIMEXP+7),
C                       CONTAINING THE PART OF THE EPSILON TABLE
C                       WICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
C                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED UP
C                       TO NOW, MULTIPLIED BY 1.5
C           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
C                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
C           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
C                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
C                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
C                       TRY TO DECREASE THE VALUE OF ERLARG.
C           LERR      - LOGICAL VARIABLE INDICATING WHETHER ABSERR HAS
C                       BEEN CHANGED IN MAIN DO-LOOP.
C           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
C                       IS NO LONGER ALLOWED (TRUE-VALUE)
C
C            MACHINE DEPENDENT CONSTANTS
C            ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQAGIE
       EPMACH = D1MACH(4)
C
C           TEST ON VALIDITY OF PARAMETERS
C           -----------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = 0.0D+00
      BLIST(1) = 0.1D+01
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      NEWFLG = .TRUE.
      IF(EPSABS.LE.0.0D+00.AND.EPSREL.LT.DMAX1(0.5D+02*EPMACH,0.5D-28))
     1  IER = 6
       IF(IER.EQ.6) GO TO 999
C
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
C           DETERMINE THE INTERVAL TO BE MAPPED ONTO (0,1).
C           IF INF = 2 THE INTEGRAL IS COMPUTED AS I = I1+I2, WHERE
C           I1 = INTEGRAL OF F OVER (-INFINITY,0),
C           I2 = INTEGRAL OF F OVER (0,+INFINITY).
C
      BOUN = BOUND
      IF(INF.EQ.2) BOUN = 0.0D+00
      CALL DQK15I(F,BOUN,INF,0.0D+00,0.1D+01,RESULT,ABSERR,
     1  DEFABS,RESABS)
C
C           TEST ON ACCURACY
C
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      DRES = DABS(RESULT)
      ERRBND = DMAX1(EPSABS,EPSREL*DRES)
      IF(ABSERR.LE.1.0D+02*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS).OR.
     1  ABSERR.EQ.0.0D+00) GO TO 130
C
C           INITIALIZATION
C           --------------
C
      UFLOW = D1MACH(1) * 0.2E+01
      LERR = .FALSE.
      CALL DEA(NEWFLG,RESULT,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      KTMIN = 0
      EXTRAP = .FALSE.
      NOEXT = .FALSE.
      IERRO = 0
      IROFF1 = 0
      IROFF2 = 0
      IROFF3 = 0
      KSGN = -1
      IF(DRES.GE.(0.1D+01-0.5D+02*EPMACH)*DEFABS) KSGN = 1
C
C           MAIN DO-LOOP
C           ------------
C
      DO 90 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        ERLAST = ERRMAX
        CALL DQK15I(F,BOUN,INF,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        CALL DQK15I(F,BOUN,INF,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2)GO TO 15
        IF(DABS(RLIST(MAXERR)-AREA12).GT.0.1D-04*DABS(AREA12)
     1  .OR.ERRO12.LT.0.99D+00*ERRMAX) GO TO 10
        IF(EXTRAP) IROFF2 = IROFF2+1
        IF(.NOT.EXTRAP) IROFF1 = IROFF1+1
   10   IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF3 = IROFF3+1
   15   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1+IROFF2.GE.10.OR.IROFF3.GE.20) IER = 2
        IF(IROFF2.GE.5) IERRO = 3
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT SOME POINTS OF THE INTEGRATION RANGE.
C
        IF(DMAX1(DABS(A1),DABS(B2)).LE.(0.1D+01+0.1D+03*EPMACH)*
     1  (DABS(A2)+0.1D+04*UFLOW)) IER = 4
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
        IF(ERROR2.GT.ERROR1) GO TO 20
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 30
   20   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   30   CALL DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
        IF(ERRSUM.LE.ERRBND) GO TO 115
        IF(IER.NE.0) GO TO 100
        IF(LAST.EQ.2) GO TO 80
        IF(NOEXT) GO TO 90
        ERLARG = ERLARG-ERLAST
        IF(DABS(B1-A1).GT.SMALL) ERLARG = ERLARG+ERRO12
        IF(EXTRAP) GO TO 40
C
C           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
C           SMALLEST INTERVAL.
C
        IF(DABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
        EXTRAP = .TRUE.
        NRMAX = 2
   40   IF(IERRO.EQ.3.OR.ERLARG.LE.ERTEST) GO TO 60
C
C           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
C           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER THE
C           LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
C
        ID = NRMAX
        JUPBND = LAST
        IF(LAST.GT.(2+LIMIT/2)) JUPBND = LIMIT+3-LAST
        DO 50 K = ID,JUPBND
          MAXERR = IORD(NRMAX)
          ERRMAX = ELIST(MAXERR)
          IF(DABS(BLIST(MAXERR)-ALIST(MAXERR)).GT.SMALL) GO TO 90
          NRMAX = NRMAX+1
   50   CONTINUE
C
C           PERFORM EXTRAPOLATION.
C
   60   CALL DEA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
        KTMIN = KTMIN+1
        IF((KTMIN.GT.5).AND.(ABSERR.LT.0.1E-02*ERRSUM).AND.(LERR))
     1     IER = 5
        IF((ABSEPS.GE.ABSERR).AND.(LERR)) GO TO 70
        KTMIN = 0
        ABSERR = ABSEPS
        LERR = .TRUE.
        RESULT = RESEPS
        CORREC = ERLARG
        ERTEST = DMAX1(EPSABS,EPSREL*DABS(RESEPS))
        IF((ABSERR.LE.ERTEST).AND.(LERR)) GO TO 100
C
C            PREPARE BISECTION OF THE SMALLEST INTERVAL.
C
   70   IF(RLIST2(LIMEXP+3).EQ.1) NOEXT = .TRUE.
        IF(IER.EQ.5) GO TO 100
        MAXERR = IORD(1)
        ERRMAX = ELIST(MAXERR)
        NRMAX = 1
        EXTRAP = .FALSE.
        SMALL = SMALL*0.5D+00
        ERLARG = ERRSUM
        GO TO 90
   80   SMALL = 0.375D+00
        ERLARG = ERRSUM
        ERTEST = ERRBND
        CALL DEA(NEWFLG,AREA,LIMEXP,RESEPS,ABSEPS,RLIST2,IERR)
   90 CONTINUE
C
C           SET FINAL RESULT AND ERROR ESTIMATE.
C           ------------------------------------
C
  100 IF(.NOT.LERR) GO TO 115
      IF((IER+IERRO).EQ.0) GO TO 110
      IF(IERRO.EQ.3) ABSERR = ABSERR+CORREC
      IF(IER.EQ.0) IER = 3
      IF(RESULT.NE.0.0D+00.AND.AREA.NE.0.0D+00)GO TO 105
      IF(ABSERR.GT.ERRSUM)GO TO 115
      IF(AREA.EQ.0.0D+00) GO TO 130
      GO TO 110
  105 IF(ABSERR/DABS(RESULT).GT.ERRSUM/DABS(AREA))GO TO 115
C
C           TEST ON DIVERGENCE
C
  110 IF(KSGN.EQ.(-1).AND.DMAX1(DABS(RESULT),DABS(AREA)).LE.
     1 DEFABS*0.1D-01) GO TO 130
      IF(0.1D-01.GT.(RESULT/AREA).OR.(RESULT/AREA).GT.0.1D+03.
     1OR.ERRSUM.GT.DABS(AREA)) IER = 6
      GO TO 130
C
C           COMPUTE GLOBAL INTEGRAL SUM.
C
  115 RESULT = 0.0D+00
      DO 120 K = 1,LAST
        RESULT = RESULT+RLIST(K)
  120 CONTINUE
      ABSERR = ERRSUM
  130 NEVAL = 30*LAST-15
      IF(INF.EQ.2) NEVAL = 2*NEVAL
      IF(IER.GT.2) IER=IER-1
  999 RETURN
      END
      SUBROUTINE DQFORM(M,N,Q,LDQ,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQFORM
C***REFER TO  DNSQ,DNSQE
C
C     SUBROUTINE DQFORM
C
C     This subroutine proceeds from the computed QR factorization of
C     an M by N matrix A to accumulate the M by M orthogonal matrix
C     Q from its factored form.
C
C     The subroutine statement is
C
C       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A and the order of Q.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       Q is an M by M array. On input the full lower trapezoid in
C         the first MIN(M,N) columns of Q contains the factored form.
C         On output Q has been accumulated into a square matrix.
C
C       LDQ is a positive integer input variable not less than M
C         which specifies the leading dimension of the array Q.
C
C       WA is a work array of length M.
C
C     Subprograms called
C
C       FORTRAN-Supplied ... MIN0
C
C     Argonne National Laboratory. MINPACK Project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DQFORM
      INTEGER MIN0
      INTEGER I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
      DOUBLE PRECISION ONE, Q(LDQ,M), SUM, TEMP, WA(M), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
C
C***FIRST EXECUTABLE STATEMENT  DQFORM
      MINMN = MIN0(M,N)
      IF (MINMN .LT. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
C
      NP1 = N + 1
      IF (M .LT. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
C
C     ACCUMULATE Q FROM ITS FACTORED FORM.
C
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .EQ. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
      END
      SUBROUTINE DQK15(F,A,B,RESULT,ABSERR,RESABS,RESASC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQK15
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  870530   (YYMMDD)
C***CATEGORY NO.  H2A1A2
C***KEYWORDS  15-POINT GAUSS-KRONROD RULES
C***AUTHOR  PIESSENS, ROBERT, AND DE DONCKER, ELISE,
C             APPLIED MATH. AND PROGR. DIV. - K. U. LEUVEN
C***PURPOSE  To compute I = Integral of F over (A,B), with error estimate,
C               and     J = integral of ABS(F) over (A,B)
C***DESCRIPTION
C           Double precision version
C
C           PARAMETERS ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              A      - Double precision: Lower limit of integration
C
C              B      - Double precision: Upper limit of integration
C
C            PARAMETERS ON RETURN
C              RESULT - Double precision: Approximation to the integral I
C                       Result is computed by applying the 15-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the7-POINT GAUSS RULE(RESG).
C
C              ABSERR - Double precision: Estimate of the modulus of the
C                       absolute error, which should not exceed ABS(I-RESULT)
C
C              RESABS - Double precision: Approximation to the integral J
C
C              RESASC - Double precision: Approximation to the integral of
C                       ABS(F-I/(B-A))  over (A,B)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DQK15
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC,
     2  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONRON QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.1294849661 6886969327 0611432679 082 D0 /
      DATA WG  (  2) / 0.2797053914 8927666790 1467771423 780 D0 /
      DATA WG  (  3) / 0.3818300505 0511894495 0369775488 975 D0 /
      DATA WG  (  4) / 0.4179591836 7346938775 5102040816 327 D0 /
C
      DATA XGK (  1) / 0.9914553711 2081263920 6854697526 329 D0 /
      DATA XGK (  2) / 0.9491079123 4275852452 6189684047 851 D0 /
      DATA XGK (  3) / 0.8648644233 5976907278 9712788640 926 D0 /
      DATA XGK (  4) / 0.7415311855 9939443986 3864773280 788 D0 /
      DATA XGK (  5) / 0.5860872354 6769113029 4144838258 730 D0 /
      DATA XGK (  6) / 0.4058451513 7739716690 6606412076 961 D0 /
      DATA XGK (  7) / 0.2077849550 0789846760 0689403773 245 D0 /
      DATA XGK (  8) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0229353220 1052922496 3732008058 970 D0 /
      DATA WGK (  2) / 0.0630920926 2997855329 0700663189 204 D0 /
      DATA WGK (  3) / 0.1047900103 2225018383 9876322541 518 D0 /
      DATA WGK (  4) / 0.1406532597 1552591874 5189590510 238 D0 /
      DATA WGK (  5) / 0.1690047266 3926790282 6583426598 550 D0 /
      DATA WGK (  6) / 0.1903505780 6478540991 3256402421 014 D0 /
      DATA WGK (  7) / 0.2044329400 7529889241 4161999234 649 D0 /
      DATA WGK (  8) / 0.2094821410 8472782801 2999174891 714 D0 /
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK15
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(8)*DABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*DMIN1(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = DMAX1
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
      SUBROUTINE DQK15I(F,BOUN,INF,A,B,RESULT,ABSERR,RESABS,RESASC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQK15I
C***DATE WRITTEN   800101   (YYMMDD)
C***REVISION DATE  830518   (YYMMDD)
C***CATEGORY NO.  H2A3A2,H2A4A2
C***KEYWORDS  15-POINT TRANSFORMED GAUSS-KRONROD RULES
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  The original (infinite integration range is mapped
C            onto the interval (0,1) and (A,B) is a part of (0,1).
C            it is the purpose to compute
C            I = Integral of transformed integrand over (A,B),
C            J = Integral of ABS(Transformed Integrand) over (A,B).
C***DESCRIPTION
C
C           Integration Rule
C           Standard Fortran subroutine
C           Double precision version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Double precision
C                       Fuction subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              BOUN   - Double precision
C                       Finite bound of original integration
C                       Range (SET TO ZERO IF INF = +2)
C
C              INF    - Integer
C                       If INF = -1, the original interval is
C                                   (-INFINITY,BOUND),
C                       If INF = +1, the original interval is
C                                   (BOUND,+INFINITY),
C                       If INF = +2, the original interval is
C                                   (-INFINITY,+INFINITY) AND
C                       The integral is computed as the sum of two
C                       integrals, one over (-INFINITY,0) and one over
C                       (0,+INFINITY).
C
C              A      - Double precision
C                       Lower limit for integration over subrange
C                       of (0,1)
C
C              B      - Double precision
C                       Upper limit for integration over subrange
C                       of (0,1)
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       Result is computed by applying the 15-POINT
C                       KRONROD RULE(RESK) obtained by optimal addition
C                       of abscissae to the 7-POINT GAUSS RULE(RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral J
C
C              RESASC - Double precision
C                       Approximation to the integral of
C                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***END PROLOGUE  DQK15I
C
      DOUBLE PRECISION A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,DABS,DINF,
     1  DMAX1,DMIN1,D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,
     2  RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,UFLOW,WG,WGK,
     3  XGK
      INTEGER INF,J
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
C           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
C           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
C                    TO THE ABSCISSAE XGK(2), XGK(4), ...
C                    WG(1), WG(3), ... ARE SET TO ZERO.
C
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
     1     0.9914553711208126D+00,     0.9491079123427585D+00,
     2     0.8648644233597691D+00,     0.7415311855993944D+00,
     3     0.5860872354676911D+00,     0.4058451513773972D+00,
     4     0.2077849550078985D+00,     0.0000000000000000D+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
     1     0.2293532201052922D-01,     0.6309209262997855D-01,
     2     0.1047900103222502D+00,     0.1406532597155259D+00,
     3     0.1690047266392679D+00,     0.1903505780647854D+00,
     4     0.2044329400752989D+00,     0.2094821410847278D+00/
C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.0000000000000000D+00,     0.1294849661688697D+00,
     2     0.0000000000000000D+00,     0.2797053914892767D+00,
     3     0.0000000000000000D+00,     0.3818300505051189D+00,
     4     0.0000000000000000D+00,     0.4179591836734694D+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC*  - ABSCISSA
C           TABSC* - TRANSFORMED ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
C                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK15I
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
      DINF = MIN0(1,INF)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      TABSC1 = BOUN+DINF*(0.1D+01-CENTR)/CENTR
      FVAL1 = F(TABSC1)
      IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
      FC = (FVAL1/CENTR)/CENTR
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ERROR.
C
      RESG = WG(8)*FC
      RESK = WGK(8)*FC
      RESABS = DABS(RESK)
      DO 10 J=1,7
        ABSC = HLGTH*XGK(J)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        TABSC1 = BOUN+DINF*(0.1D+01-ABSC1)/ABSC1
        TABSC2 = BOUN+DINF*(0.1D+01-ABSC2)/ABSC2
        FVAL1 = F(TABSC1)
        FVAL2 = F(TABSC2)
        IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
        IF(INF.EQ.2) FVAL2 = FVAL2+F(-TABSC2)
        FVAL1 = (FVAL1/ABSC1)/ABSC1
        FVAL2 = (FVAL2/ABSC2)/ABSC2
        FV1(J) = FVAL1
        FV2(J) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(J)*FSUM
        RESABS = RESABS+WGK(J)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(8)*DABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESASC = RESASC*HLGTH
      RESABS = RESABS*HLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.D0) ABSERR = RESASC*
     1 DMIN1(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = DMAX1
     1 ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
      SUBROUTINE DQPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQPSRT
C***REFER TO  DQAGE,DQAGIE,DQAGPE,DQAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION DATE  810101   (YYMMDD)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  PIESSENS, ROBERT, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C           DE DONCKER, ELISE, APPLIED MATH. AND PROGR. DIV. -
C             K. U. LEUVEN
C***PURPOSE  This routine maintains the descending ordering in the
C            list of the local error estimated resulting from the
C            interval subdivision process. At each call two error
C            estimates are inserted using the sequential search
C            method, top-down for the largest error estimate and
C            bottom-up for the smallest error estimate.
C***DESCRIPTION
C
C           Ordering routine
C           Standard fortran subroutine
C           Double precision version
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - Integer
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - Integer
C                       Number of error estimates currently in the list
C
C              MAXERR - Integer
C                       Maxerr points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - Double precision
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - Double precision
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - Integer
C                       Vector of dimension LAST, the first K elements
C                       of which contain pointers to the error
C                       estimates, such that
C                       ELIST(IORD(1)),...,  ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - Integer
C                       MAXERR = IORD(NRMAX)
C***END PROLOGUE  DQPSRT
C
      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(LAST),IORD(LAST)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  DQPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
C           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
C           SUBDIVISIONS STILL ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
      SUBROUTINE DQRANK(A,LDA,M,N,TOL,KR,JPVT,QRAUX,WORK)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQRANK
C***REVISION NOVEMBER 15, 1980
C***CATEGORY NO.  E2G1A
C***KEYWORD(S)  OVERDETERMINED,LEAST SQUARES,LINEAR EQUATIONS
C***AUTHOR  D. KAHANER (NBS)
C***DATE WRITTEN
C***PURPOSE
C      COMPUTES THE QR FACTORIZATION OF AN  M  BY  N  MATRIX A.  THIS 
C      ROUTINE IS USED IN CONJUNCTION WITH DQRLSS TO SOLVE LINEAR
C      SYSTEMS OF EQUATIONS IN A LEAST SQUARE SENSE.
C***DESCRIPTION
C
C     DQRANK IS USED IN CONJUNCTION WITH DQRLSS TO SOLVE
C     OVERDETERMINED, UNDERDETERMINED AND SINGULAR LINEAR SYSTEMS
C     IN A LEAST SQUARES SENSE.
C     DQRANK USES THE LINPACK SUBROUTINE DQRDC TO COMPUTE THE QR
C     FACTORIZATION, WITH COLUMN PIVOTING, OF AN  M  BY  N  MATRIX  A .
C     THE NUMERICAL RANK IS DETERMINED USING THE TOLERANCE TOL.
C
C     FOR MORE INFORMATION, SEE CHAPTER 9 OF LINPACK USERS GUIDE,
C     J. DONGARRA ET ALL, SIAM, 1979.
C
C     ON ENTRY
C
C        A     DOUBLE PRECISION (LDA,N) . 
C              THE MATRIX WHOSE DECOMPOSITION IS TO BE COMPUTED.
C
C        LDA   INTEGER.
C              THE LEADING DIMENSION OF A .
C
C        M     INTEGER.
C              THE NUMBER OF ROWS OF A .
C
C        N     INTEGER.
C              THE NUMBER OF COLUMNS OF  A .
C
C        TOL   DOUBLE PRECISION.
C              A RELATIVE TOLERANCE USED TO DETERMINE THE NUMERICAL
C              RANK.  THE PROBLEM SHOULD BE SCALED SO THAT ALL THE 
C              ELEMENTS OF  A   HAVE ROUGHLY THE SAME ABSOLUTE ACCURACY
C              EPS.  THEN A REASONABLE VALUE FOR  TOL  IS ROUGHLY  EPS  
C              DIVIDED BY THE MAGNITUDE OF THE LARGEST ELEMENT.
C
C        JPVT  INTEGER(N)
C
C        QRAUX DOUBLE PRECISION(N)
C
C        WORK  DOUBLE PRECISION(N)
C              THREE AUXILLIARY VECTORS USED BY DQRDC .
C
C     ON RETURN
C
C        A     CONTAINS THE OUTPUT FROM DQRDC.
C              THE TRIANGULAR MATRIX  R  OF THE QR FACTORIZATION IS
C              CONTAINED IN THE UPPER TRIANGLE AND INFORMATION NEEDED 
C              TO RECOVER THE ORTHOGONAL MATRIX  Q  IS STORED BELOW
C              THE DIAGONAL IN  A  AND IN THE VECTOR QRAUX .
C
C        KR    INTEGER.
C              THE NUMERICAL RANK.
C
C        JPVT  THE PIVOT INFORMATION FROM DQRDC.
C
C     COLUMNS JPVT(1),...,JPVT(KR) OF THE ORIGINAL MATRIX ARE LINEARLY
C     INDEPENDENT TO WITHIN THE TOLERANCE TOL AND THE REMAINING COLUMNS
C     ARE LINEARLY DEPENDENT.  ABS(A(1,1))/ABS(A(KR,KR))  IS AN ESTIMATE
C     OF THE CONDITION NUMBER OF THE MATRIX OF INDEPENDENT COLUMNS,
C     AND OF  R .  THIS ESTIMATE WILL BE .LE. 1/TOL .
C
C      USAGE.....SEE SUBROUTINE DQRLSS
C
C***REFERENCE(S)
C      DONGARRA, ET AL, LINPACK USERS GUIDE, SIAM, 1979
C***ROUTINES CALLED  DQRDC
C***END PROLOGUE
      INTEGER LDA,M,N,KR,JPVT(1)
      DOUBLE PRECISION A(LDA,1),TOL,QRAUX(1),WORK(1)
      INTEGER J,K
C***FIRST EXECUTABLE STATEMENT
      DO 10 J = 1, N
         JPVT(J) = 0
   10 CONTINUE
      CALL DQRDC(A,LDA,M,N,QRAUX,JPVT,WORK,1)
      KR = 0
      K = MIN0(M,N) 
      DO 20 J = 1, K
         IF (ABS(A(J,J)) .LE. TOL*ABS(A(1,1))) GO TO 30
         KR = J
   20 CONTINUE
   30 RETURN
      END 
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)

c*********************************************************************72
c            
C***BEGIN PROLOGUE  DQRDC                                               
C***DATE WRITTEN   780814   (YYMMDD)                                    
C***REVISION DATE  861211   (YYMMDD)                                    
C***CATEGORY NO.  D5                                                    
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),                                  
C             TYPE=DOUBLE PRECISION(SQRDC-S DQRDC-D CQRDC-C),           
C             LINEAR ALGEBRA,MATRIX,ORTHOGONAL TRIANGULAR,              
C             QR DECOMPOSITION                                          
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)                            
C***PURPOSE  Uses Householder transformations to compute the Qr factori-
C            zation of N by P matrix X.  Column pivoting is optional.   
C***DESCRIPTION                                                         
C                                                                       
C     DQRDC uses Householder transformations to compute the QR          
C     factorization of an N by P matrix X.  Column pivoting             
C     based on the 2-norms of the reduced columns may be                
C     performed at the user's option.                                   
C                                                                       
C     On Entry                                                          
C                                                                       
C        X       DOUBLE PRECISION(LDX,P), where LDX .GE. N.             
C                X contains the matrix whose decomposition is to be     
C                computed.                                              
C                                                                       
C        LDX     INTEGER.                                               
C                LDX is the leading dimension of the array X.           
C                                                                       
C        N       INTEGER.                                               
C                N is the number of rows of the matrix X.               
C                                                                       
C        P       INTEGER.                                               
C                P is the number of columns of the matrix X.            
C                                                                       
C        JPVT    INTEGER(P).                                            
C                JPVT contains integers that control the selection      
C                of the pivot columns.  The K-th column X(K) of X       
C                is placed in one of three classes according to the     
C                value of JPVT(K).                                      
C                                                                       
C                   If JPVT(K) .GT. 0, then X(K) is an initial          
C                                      column.                          
C                                                                       
C                   If JPVT(K) .EQ. 0, then X(K) is a free column.      
C                                                                       
C                   If JPVT(K) .LT. 0, then X(K) is a final column.     
C                                                                       
C                Before the decomposition is computed, initial columns  
C                are moved to the beginning of the array X and final    
C                columns to the end.  Both initial and final columns    
C                are frozen in place during the computation and only    
C                free columns are moved.  At the K-th stage of the      
C                reduction, if X(K) is occupied by a free column        
C                it is interchanged with the free column of largest     
C                reduced norm.  JPVT is not referenced if               
C                JOB .EQ. 0.                                            
C                                                                       
C        WORK    DOUBLE PRECISION(P).                                   
C                WORK is a work array.  WORK is not referenced if       
C                JOB .EQ. 0.                                            
C                                                                       
C        JOB     INTEGER.                                               
C                JOB is an integer that initiates column pivoting.      
C                If JOB .EQ. 0, no pivoting is done.                    
C                If JOB .NE. 0, pivoting is done.                       
C                                                                       
C     On Return                                                         
C                                                                       
C        X       X contains in its upper triangle the upper             
C                triangular matrix R of the QR factorization.           
C                Below its diagonal X contains information from         
C                which the orthogonal part of the decomposition         
C                can be recovered.  Note that if pivoting has           
C                been requested, the decomposition is not that          
C                of the original matrix X but that of X                 
C                with its columns permuted as described by JPVT.        
C                                                                       
C        QRAUX   DOUBLE PRECISION(P).                                   
C                QRAUX contains further information required to recover 
C                the orthogonal part of the decomposition.              
C                                                                       
C        JPVT    JPVT(K) contains the index of the column of the        
C                original matrix that has been interchanged into        
C                the K-th column, if pivoting was requested.            
C                                                                       
C     LINPACK.  This version dated 08/14/78 .                           
C     G. W. Stewart, University of Maryland, Argonne National Lab.      
C                                                                       
C     DQRDC uses the following functions and subprograms.               
C                                                                       
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2                                 
C     Fortran DABS,DMAX1,MIN0,DSQRT                                     
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   
C***ROUTINES CALLED  DAXPY,DDOT,DNRM2,DSCAL,DSWAP                       
C***END PROLOGUE  DQRDC                                                 
      INTEGER LDX,N,P,JOB                                               
      INTEGER JPVT(1)                                                   
      DOUBLE PRECISION X(LDX,1),QRAUX(1),WORK(1)                        
C                                                                       
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU                                 
      DOUBLE PRECISION MAXNRM,DNRM2,TT                                  
      DOUBLE PRECISION DDOT,NRMXL,T                                     
      LOGICAL NEGJ,SWAPJ                                                
C                                                                       
C***FIRST EXECUTABLE STATEMENT  DQRDC                                   
      PL = 1                                                            
      PU = 0                                                            
      IF (JOB .EQ. 0) GO TO 60                                          
C                                                                       
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS            
C        ACCORDING TO JPVT.                                             
C                                                                       
         DO 20 J = 1, P                                                 
            SWAPJ = JPVT(J) .GT. 0                                      
            NEGJ = JPVT(J) .LT. 0                                       
            JPVT(J) = J                                                 
            IF (NEGJ) JPVT(J) = -J                                      
            IF (.NOT.SWAPJ) GO TO 10                                    
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)          
               JPVT(J) = JPVT(PL)                                       
               JPVT(PL) = J                                             
               PL = PL + 1                                              
   10       CONTINUE                                                    
   20    CONTINUE                                                       
         PU = P                                                         
         DO 50 JJ = 1, P                                                
            J = P - JJ + 1                                              
            IF (JPVT(J) .GE. 0) GO TO 40                                
               JPVT(J) = -JPVT(J)                                       
               IF (J .EQ. PU) GO TO 30                                  
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)                      
                  JP = JPVT(PU)                                         
                  JPVT(PU) = JPVT(J)                                    
                  JPVT(J) = JP                                          
   30          CONTINUE                                                 
               PU = PU - 1                                              
   40       CONTINUE                                                    
   50    CONTINUE                                                       
   60 CONTINUE                                                          
C                                                                       
C     COMPUTE THE NORMS OF THE FREE COLUMNS.                            
C                                                                       
      IF (PU .LT. PL) GO TO 80                                          
      DO 70 J = PL, PU                                                  
         QRAUX(J) = DNRM2(N,X(1,J),1)                                   
         WORK(J) = QRAUX(J)                                             
   70 CONTINUE                                                          
   80 CONTINUE                                                          
C                                                                       
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.                           
C                                                                       
      LUP = MIN0(N,P)                                                   
      DO 200 L = 1, LUP                                                 
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120                        
C                                                                       
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT              
C           INTO THE PIVOT POSITION.                                    
C                                                                       
            MAXNRM = 0.0D0                                              
            MAXJ = L                                                    
            DO 100 J = L, PU                                            
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90                       
                  MAXNRM = QRAUX(J)                                     
                  MAXJ = J                                              
   90          CONTINUE                                                 
  100       CONTINUE                                                    
            IF (MAXJ .EQ. L) GO TO 110                                  
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)                       
               QRAUX(MAXJ) = QRAUX(L)                                   
               WORK(MAXJ) = WORK(L)                                     
               JP = JPVT(MAXJ)                                          
               JPVT(MAXJ) = JPVT(L)                                     
               JPVT(L) = JP                                             
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         QRAUX(L) = 0.0D0                                               
         IF (L .EQ. N) GO TO 190                                        
C                                                                       
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.        
C                                                                       
            NRMXL = DNRM2(N-L+1,X(L,L),1)                               
            IF (NRMXL .EQ. 0.0D0) GO TO 180                             
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))       
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)                   
               X(L,L) = 1.0D0 + X(L,L)                                  
C                                                                       
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,       
C              UPDATING THE NORMS.                                      
C                                                                       
               LP1 = L + 1                                              
               IF (P .LT. LP1) GO TO 170                                
               DO 160 J = LP1, P                                        
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)             
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)                 
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150               
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150                    
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2            
                     TT = DMAX1(TT,0.0D0)                               
                     T = TT                                             
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2       
                     IF (TT .EQ. 1.0D0) GO TO 130                       
                        QRAUX(J) = QRAUX(J)*DSQRT(T)                    
                     GO TO 140                                          
  130                CONTINUE                                           
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)                
                        WORK(J) = QRAUX(J)                              
  140                CONTINUE                                           
  150             CONTINUE                                              
  160          CONTINUE                                                 
  170          CONTINUE                                                 
C                                                                       
C              SAVE THE TRANSFORMATION.                                 
C                                                                       
               QRAUX(L) = X(L,L)                                        
               X(L,L) = -NRMXL                                          
  180       CONTINUE                                                    
  190    CONTINUE                                                       
  200 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQRFAC
C***REFER TO  DNLS1,DNLS1E,DNSQ,DNSQE
C***ROUTINES CALLED  D1MACH,DENORM
C***DESCRIPTION
C
C   **** Double Precision version of QRFAC ****
C
C     **********
C     Subroutine DQRFAC
C
C     This subroutine uses Householder transformations with column
C     pivoting (optional) to compute a QR factorization of the
C     M by N matrix A. That is, DQRFAC determines an orthogonal
C     matrix Q, a permutation matrix P, and an upper trapezoidal
C     matrix R with diagonal elements of nonincreasing magnitude,
C     such that A*P = Q*R. The Householder transformation for
C     column K, K = 1,2,...,MIN(M,N), is of the form
C
C                           T
C           I - (1/U(K))*U*U
C
C     where U has zeros in the first K-1 positions. The form of
C     this transformation and the method of pivoting first
C     appeared in the corresponding LINPACK subroutine.
C
C     The subroutine statement is
C
C       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       A is an M by N array. On input A contains the matrix for
C         which the QR factorization is to be computed. On output
C         the strict upper trapezoidal part of A contains the strict
C         upper trapezoidal part of R, and the lower trapezoidal
C         part of A contains a factored form of Q (the non-trivial
C         elements of the U vectors described above).
C
C       LDA is a positive integer input variable not less than M
C         which specifies the leading dimension of the array A.
C
C       PIVOT is a logical input variable. If pivot is set .TRUE.,
C         then column pivoting is enforced. If pivot is set .FALSE.,
C         then no column pivoting is done.
C
C       IPVT is an integer output array of length LIPVT. IPVT
C         defines the permutation matrix P such that A*P = Q*R.
C         Column J of P is column IPVT(J) of the identity matrix.
C         If pivot is .FALSE., IPVT is not referenced.
C
C       LIPVT is a positive integer input variable. If PIVOT is
C             .FALSE., then LIPVT may be as small as 1. If PIVOT is
C             .TRUE., then LIPVT must be at least N.
C
C       SIGMA is an output array of length N which contains the
C         diagonal elements of R.
C
C       ACNORM is an output array of length N which contains the
C         norms of the corresponding columns of the input matrix A.
C         If this information is not needed, then ACNORM can coincide
C         with SIGMA.
C
C       WA is a work array of length N. If pivot is .FALSE., then WA
C         can coincide with SIGMA.
C
C     Subprograms called
C
C       SLATEC-Supplied ... D1MACH,DENORM
C       FORTRAN-Supplied ... DMAX1,DSQRT,MIN0
C
C     MINPACK. Version of December 1978.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C     **********
C***END PROLOGUE  DQRFAC
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      SAVE ONE, P05, ZERO
      DOUBLE PRECISION A(LDA,N),SIGMA(N),ACNORM(N),WA(N)
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION D1MACH,DENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
C***FIRST EXECUTABLE STATEMENT  DQRFAC
      EPSMCH = D1MACH(4)
C
C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
         ACNORM(J) = DENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
         AJNORM = DENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
C
C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C        AND UPDATE THE NORMS.
C
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            SIGMA(K) = DENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DQRFAC.
C
      END
      SUBROUTINE DQRLS (A,LDA,M,N,TOL,KR,B,X,RSD,JPVT,QRAUX,WORK,
     *                  ITASK,IND)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQRLS
C***REVISION SEPTEMBER 11, 1987
C***CATEGORY NO.  E2G1A
C***KEYWORD(S)  OVERDETERMINED,LEAST SQUARES,LINEAR EQUATIONS
C***AUTHOR  STEPHEN NASH (GEORGE MASON UNIVERSITY)
C***DATE WRITTEN SEPTEMBER 11, 1987
C***PURPOSE
C      SOLVES AN OVERDETERMINED, UNDERDETERMINED, OR SINGULAR SYSTEM OF
C      LINEAR EQUATIONS IN LEAST SQUARE SENSE.  THE SOLUTION IS OBTAINED
C      USING A QR FACTORIZATION OF THE  M  BY  N  COEFFICIENT MATRIX  A.
C
C***DESCRIPTION
C     DQRLS IS USED TO SOLVE OVERDETERMINED, UNDERDETERMINED AND 
C     SINGULAR LINEAR SYSTEMS IN A LEAST SQUARES SENSE.
C     THE SYSTEM IS  A*X  APPROXIMATES  B  WHERE  A  IS  M  BY  N.
C     B  IS A GIVEN  M-VECTOR, AND  X  IS THE  N-VECTOR TO BE COMPUTED.  
C     A SOLUTION  X  IS FOUND WHICH MINIMIZES THE SUM OF SQUARES (2-NORM)
C     OF THE RESIDUAL,  A*X - B .
C     DQRLS USES THE LINPACK SUBROUTINE SQRDC TO COMPUTE THE QR
C     FACTORIZATION, WITH COLUMN PIVOTING, OF AN  M  BY  N  MATRIX  A .
C     THE NUMERICAL RANK IS DETERMINED USING THE TOLERANCE TOL.
C
C     FOR MORE INFORMATION, SEE CHAPTER 9 OF LINPACK USERS GUIDE,
C     J. DONGARRA ET ALL, SIAM, 1979.
C
C
C     ON ENTRY
C
C        A     DOUBLE PRECISION (LDA,N) . 
C              THE MATRIX WHOSE DECOMPOSITION IS TO BE COMPUTED.
C
C        LDA   INTEGER.
C              THE LEADING DIMENSION OF A .
C              M <= LDA.
C
C        M     INTEGER.
C              THE NUMBER OF ROWS OF A .
C
C        N     INTEGER.
C              THE NUMBER OF COLUMNS OF  A .
C
C        TOL   DOUBLE PRECISION.
C              A RELATIVE TOLERANCE USED TO DETERMINE THE NUMERICAL
C              RANK.  THE PROBLEM SHOULD BE SCALED SO THAT ALL THE 
C              ELEMENTS OF  A   HAVE ROUGHLY THE SAME ABSOLUTE ACCURACY
C              EPS.  THEN A REASONABLE VALUE FOR  TOL  IS ROUGHLY  EPS  
C              DIVIDED BY THE MAGNITUDE OF THE LARGEST ELEMENT.
C
C        JPVT  INTEGER(N)
C
C        QRAUX DOUBLE PRECISION(N)
C
C        WORK  DOUBLE PRECISION(N)
C              THREE AUXILIARY VECTORS USED TO FACTOR THE MATRIX A.
C              (ARRAY WORK IS NOT REQUIRED IF ITASK .GT. 1)
C
C        B     DOUBLE PRECISION(M)
C              THE RIGHT HAND SIDE OF THE LINEAR SYSTEM.
C
C        ITASK INTEGER.
C              IF ITASK=1, THEN DQRLS FACTORS THE MATRIX A AND
C                          SOLVES THE LEAST SQUARES PROBLEM.
C              IF ITASK>1, THEN DQRLS ASSUMES THAT THE MATRIX A
C                          WAS FACTORED WITH AN EARLIER CALL TO
C                          DQRLS, AND ONLY SOLVES THE LEAST SQUARES
C                          PROBLEM.
C
C     ON RETURN
C
C        X     DOUBLE PRECISION(N) .
C              A LEAST SQUARES SOLUTION TO THE LINEAR SYSTEM.
C
C        RSD   DOUBLE PRECISION(M) .
C              THE RESIDUAL, B - A*X .  RSD MAY OVERWRITE  B .
C
C        IND   INTEGER
C              ERROR CODE:  IND =  0:  NO ERROR
C                           IND = -1: N .GT. LDA   (FATAL ERROR)
C                           IND = -2: N .LT. 1     (FATAL ERROR)
C                           IND = -3: ITASK .LT. 1 (FATAL ERROR)
C
C        A     CONTAINS THE OUTPUT FROM SQRDC.
C              THE TRIANGULAR MATRIX  R  OF THE QR FACTORIZATION IS
C              CONTAINED IN THE UPPER TRIANGLE AND INFORMATION NEEDED 
C              TO RECOVER THE ORTHOGONAL MATRIX  Q  IS STORED BELOW
C              THE DIAGONAL IN  A  AND IN THE VECTOR QRAUX .
C
C        KR    INTEGER.
C              THE NUMERICAL RANK.
C
C        JPVT  THE PIVOT INFORMATION FROM SQRDC.
C
C     COLUMNS JPVT(1),...,JPVT(KR) OF THE ORIGINAL MATRIX ARE LINEARLY
C     INDEPENDENT TO WITHIN THE TOLERANCE TOL AND THE REMAINING COLUMNS
C     ARE LINEARLY DEPENDENT.  ABS(A(1,1))/ABS(A(KR,KR))  IS AN ESTIMATE
C     OF THE CONDITION NUMBER OF THE MATRIX OF INDEPENDENT COLUMNS,
C     AND OF  R .  THIS ESTIMATE WILL BE .LE. 1/TOL .
C
C      USAGE....
C        DQRLS CAN BE EFFICIENTLY USED TO SOLVE SEVERAL LEAST SQUARES
C      PROBLEMS WITH THE SAME MATRIX A.  THE FIRST SYSTEM IS SOLVED
C      WITH ITASK = 1.  THE SUBSEQUENT SYSTEMS ARE SOLVED WITH
C      ITASK = 2, TO AVOID THE RECOMPUTATION OF THE MATRIX FACTORS.
C      THE PARAMETERS  KR, JPVT, AND QRAUX MUST NOT BE MODIFIED
C      BETWEEN CALLS TO DQRLS.
C
C***REFERENCE(S)
C      DONGARRA, ET AL, LINPACK USERS GUIDE, SIAM, 1979
C***ROUTINES CALLED  DQRANK, DQRLSS, XERROR
C***END PROLOGUE
      INTEGER  LDA, M, N, ITASK, JPVT(1), KR, IND
      DOUBLE PRECISION  A(LDA,1), B(1), X(1), RSD(1), QRAUX(1), 
     *     WORK(1), TOL
      CHARACTER  MSG*54
C***FIRST EXECUTABLE STATEMENT
      IF (LDA   .LT. M) GO TO 10
      IF (N     .LE. 0) GO TO 20
      IF (ITASK .LT. 1) GO TO 30
      IND = 0
C
      IF (ITASK .EQ. 1)
C
C FACTOR MATRIX
C
     *    CALL DQRANK (A, LDA, M, N, TOL, KR, JPVT, QRAUX, WORK)
C
C SOLVE LEAST-SQUARES PROBLEM
C
      CALL DQRLSS (A, LDA, M, N, KR, B, X, RSD, JPVT, QRAUX)
      RETURN
C
C ERROR IN CALLING SEQUENCE
C
C LDA .LT. N
10    IND = -1
      WRITE (MSG, '(
     *  ''DQRLS ERROR (IND=-1) -- LDA='', I5, '' IS LESS THAN M='',
     *      I5       )'   ) LDA, M
      CALL XERROR(MSG(1:54), 54, -1, 0)
      RETURN
C
C N .LT. 1
20    IND = -2
      WRITE (MSG, '(
     *  ''DQRLS ERROR (IND=-2) -- N='', I5, '' IS LESS THAN 1.'') ')N
      CALL XERROR(MSG(1:47), 47, -2, 0)
      RETURN
C
C ITASK .LT. 1
30    IND = -3
      WRITE (MSG, '(
     *  ''DQRLS ERROR (IND=-3) -- ITASK='', I5, '' IS LESS THAN 1.'')
     *               ') ITASK
      CALL XERROR(MSG(1:51), 51, -3, 0)
      RETURN
      END
      SUBROUTINE DQRLSS(A,LDA,M,N,KR,B,X,RSD,JPVT,QRAUX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DQRLSS
C***REVISION NOVEMBER 15, 1980
C***CATEGORY NO.  E2G1A
C***KEYWORD(S)  OVERDETERMINED,LEAST SQUARES,LINEAR EQUATIONS
C***AUTHOR  D. KAHANER (NBS)
C***DATE WRITTEN
C***PURPOSE
C      SOLVES AN OVERDETERMINED, UNDERDETERMINED, OR SINGULAR SYSTEM OF
C      LINEAR EQUATIONS IN LEAST SQUARE SENSE.
C***DESCRIPTION
C
C     DQRLSS USES THE LINPACK SUBROUTINE DQRSL TO SOLVE AN 
C     OVERDETERMINED, UNDERDETERMINED, OR SINGULAR LINEAR SYSTEM IN A 
C     LEAST SQUARES SENSE.  IT MUST BE PRECEEDED BY DQRANK .
C     THE SYSTEM IS  A*X  APPROXIMATES  B  WHERE  A  IS  M  BY  N  WITH
C     RANK  KR  (DETERMINED BY DQRANK),  B  IS A GIVEN  M-VECTOR,
C     AND  X  IS THE  N-VECTOR TO BE COMPUTED.  A SOLUTION  X  WITH
C     AT MOST  KR  NONZERO COMPONENTS IS FOUND WHICH MINIMIMZES
C     THE 2-NORM OF THE RESIDUAL,  A*X - B .
C
C     ON ENTRY
C
C        A,LDA,M,N,KR,JPVT,QRAUX
C              THE OUTPUT FROM DQRANK . 
C
C        B     DOUBLE PRECISION(M) .
C              THE RIGHT HAND SIDE OF THE LINEAR SYSTEM.
C
C     ON RETURN
C
C        X     DOUBLE PRECISION(N) .
C              A LEAST SQUARES SOLUTION TO THE LINEAR SYSTEM.
C
C        RSD   DOUBLE PRECISION(M) .
C              THE RESIDUAL, B - A*X .  RSD MAY OVERWITE  B .
C
C      USAGE....
C        ONCE THE MATRIX A HAS BEEN FORMED, DQRANK SHOULD BE
C      CALLED ONCE TO DECOMPOSE IT.  THEN FOR EACH RIGHT HAND
C      SIDE, B, DQRLSS SHOULD BE CALLED ONCE TO OBTAIN THE
C      SOLUTION AND RESIDUAL. 
C
C     INTEGER J,K,INFO
C     DOUBLE PRECISION T
C
C***REFERENCE(S)
C      DONGARRA, ET AL, LINPACK USERS GUIDE, SIAM, 1979
C***ROUTINES CALLED  DQRSL
C***END PROLOGUE
      INTEGER LDA,M,N,KR,JPVT(1)
      DOUBLE PRECISION A(LDA,1),B(1),X(1),RSD(1),QRAUX(1)
C***FIRST EXECUTABLE STATEMENT
      IF (KR .NE. 0)
     1   CALL DQRSL(A,LDA,M,KR,QRAUX,B,RSD,RSD,X,RSD,RSD,110,INFO)
      DO 40 J = 1, N
         JPVT(J) = -JPVT(J)
         IF (J .GT. KR) X(J) = 0.0D0
   40 CONTINUE
      DO 70 J = 1, N
         IF (JPVT(J) .GT. 0) GO TO 70
         K = -JPVT(J)
         JPVT(J) = K
   50    CONTINUE
            IF (K .EQ. J) GO TO 60
            T = X(J)
            X(J) = X(K)
            X(K) = T
            JPVT(K) = -JPVT(K)
            K = JPVT(K)
         GO TO 50
   60    CONTINUE
   70 CONTINUE
      RETURN
      END 
      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)   

c*********************************************************************72
c   
C***BEGIN PROLOGUE  DQRSL                                               
C***DATE WRITTEN   780814   (YYMMDD)                                    
C***REVISION DATE  861211   (YYMMDD)                                    
C***CATEGORY NO.  D9,D2A1                                               
C***KEYWORDS  LIBRARY=SLATEC(LINPACK),                                  
C             TYPE=DOUBLE PRECISION(SQRSL-S DQRSL-D CQRSL-C),           
C             LINEAR ALGEBRA,MATRIX,ORTHOGONAL TRIANGULAR,SOLVE         
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)                            
C***PURPOSE  Applies the output of DQRDC to compute coordinate          
C            transformations, projections, and least squares solutions. 
C***DESCRIPTION                                                         
C                                                                       
C     DQRSL applies the output of DQRDC to compute coordinate           
C     transformations, projections, and least squares solutions.        
C     For K .LE. MIN(N,P), let XK be the matrix                         
C                                                                       
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))              
C                                                                       
C     formed from columnns JPVT(1), ... ,JPVT(K) of the original        
C     N X P matrix X that was input to DQRDC (if no pivoting was        
C     done, XK consists of the first K columns of X in their            
C     original order).  DQRDC produces a factored orthogonal matrix Q   
C     and an upper triangular matrix R such that                        
C                                                                       
C              XK = Q * (R)                                             
C                       (0)                                             
C                                                                       
C     This information is contained in coded form in the arrays         
C     X and QRAUX.                                                      
C                                                                       
C     On Entry                                                          
C                                                                       
C        X      DOUBLE PRECISION(LDX,P).                                
C               X contains the output of DQRDC.                         
C                                                                       
C        LDX    INTEGER.                                                
C               LDX is the leading dimension of the array X.            
C                                                                       
C        N      INTEGER.                                                
C               N is the number of rows of the matrix XK.  It must      
C               have the same value as N in DQRDC.                      
C                                                                       
C        K      INTEGER.                                                
C               K is the number of columns of the matrix XK.  K         
C               must not be greater than MIN(N,P), where P is the       
C               same as in the calling sequence to DQRDC.               
C                                                                       
C        QRAUX  DOUBLE PRECISION(P).                                    
C               QRAUX contains the auxiliary output from DQRDC.         
C                                                                       
C        Y      DOUBLE PRECISION(N)                                     
C               Y contains an N-vector that is to be manipulated        
C               by DQRSL.                                               
C                                                                       
C        JOB    INTEGER.                                                
C               JOB specifies what is to be computed.  JOB has          
C               the decimal expansion ABCDE, with the following         
C               meaning.                                                
C                                                                       
C                    If A .NE. 0, compute QY.                           
C                    If B,C,D, or E .NE. 0, compute QTY.                
C                    If C .NE. 0, compute B.                            
C                    If D .NE. 0, compute RSD.                          
C                    If E .NE. 0, compute XB.                           
C                                                                       
C               Note that a request to compute B, RSD, or XB            
C               automatically triggers the computation of QTY, for      
C               which an array must be provided in the calling          
C               sequence.                                               
C                                                                       
C     On Return                                                         
C                                                                       
C        QY     DOUBLE PRECISION(N).                                    
C               QY contains Q*Y, if its computation has been            
C               requested.                                              
C                                                                       
C        QTY    DOUBLE PRECISION(N).                                    
C               QTY contains TRANS(Q)*Y, if its computation has         
C               been requested.  Here TRANS(Q) is the                   
C               transpose of the matrix Q.                              
C                                                                       
C        B      DOUBLE PRECISION(K)                                     
C               B contains the solution of the least squares problem    
C                                                                       
C                    minimize norm2(Y - XK*B),                          
C                                                                       
C               if its computation has been requested.  (Note that      
C               if pivoting was requested in DQRDC, the J-th            
C               component of B will be associated with column JPVT(J)   
C               of the original matrix X that was input into DQRDC.)    
C                                                                       
C        RSD    DOUBLE PRECISION(N).                                    
C               RSD contains the least squares residual Y - XK*B,       
C               if its computation has been requested.  RSD is          
C               also the orthogonal projection of Y onto the            
C               orthogonal complement of the column space of XK.        
C                                                                       
C        XB     DOUBLE PRECISION(N).                                    
C               XB contains the least squares approximation XK*B,       
C               if its computation has been requested.  XB is also      
C               the orthogonal projection of Y onto the column space    
C               of X.                                                   
C                                                                       
C        INFO   INTEGER.                                                
C               INFO is zero unless the computation of B has            
C               been requested and R is exactly singular.  In           
C               this case, INFO is the index of the first zero          
C               diagonal element of R and B is left unaltered.          
C                                                                       
C     The parameters QY, QTY, B, RSD, and XB are not referenced         
C     if their computation is not requested and in this case            
C     can be replaced by dummy variables in the calling program.        
C     To save storage, the user may in some cases use the same          
C     array for different parameters in the calling sequence.  A        
C     frequently occuring example is when one wishes to compute         
C     any of B, RSD, or XB and does not need Y or QTY.  In this         
C     case one may identify Y, QTY, and one of B, RSD, or XB, while     
C     providing separate arrays for anything else that is to be         
C     computed.  Thus the calling sequence                              
C                                                                       
C          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)         
C                                                                       
C     will result in the computation of B and RSD, with RSD             
C     overwriting Y.  More generally, each item in the following        
C     list contains groups of permissible identifications for           
C     a single calling sequence.                                        
C                                                                       
C          1. (Y,QTY,B) (RSD) (XB) (QY)                                 
C                                                                       
C          2. (Y,QTY,RSD) (B) (XB) (QY)                                 
C                                                                       
C          3. (Y,QTY,XB) (B) (RSD) (QY)                                 
C                                                                       
C          4. (Y,QY) (QTY,B) (RSD) (XB)                                 
C                                                                       
C          5. (Y,QY) (QTY,RSD) (B) (XB)                                 
C                                                                       
C          6. (Y,QY) (QTY,XB) (B) (RSD)                                 
C                                                                       
C     In any group the value returned in the array allocated to         
C     the group corresponds to the last member of the group.            
C                                                                       
C     LINPACK.  This version dated 08/14/78 .                           
C     G. W. Stewart, University of Maryland, Argonne National Lab.      
C                                                                       
C     DQRSL uses the following functions and subprograms.               
C                                                                       
C     BLAS DAXPY,DCOPY,DDOT                                             
C     Fortran DABS,MIN0,MOD                                             
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,    
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.                   
C***ROUTINES CALLED  DAXPY,DCOPY,DDOT                                   
C***END PROLOGUE  DQRSL                                                 
      INTEGER LDX,N,K,JOB,INFO                                          
      DOUBLE PRECISION X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1), 
     1                 XB(1)                                            
C                                                                       
      INTEGER I,J,JJ,JU,KP1                                             
      DOUBLE PRECISION DDOT,T,TEMP                                      
      LOGICAL CB,CQY,CQTY,CR,CXB                                        
C                                                                       
C     SET INFO FLAG.                                                    
C                                                                       
C***FIRST EXECUTABLE STATEMENT  DQRSL                                   
      INFO = 0                                                          
C                                                                       
C     DETERMINE WHAT IS TO BE COMPUTED.                                 
C                                                                       
      CQY = JOB/10000 .NE. 0                                            
      CQTY = MOD(JOB,10000) .NE. 0                                      
      CB = MOD(JOB,1000)/100 .NE. 0                                     
      CR = MOD(JOB,100)/10 .NE. 0                                       
      CXB = MOD(JOB,10) .NE. 0                                          
      JU = MIN0(K,N-1)                                                  
C                                                                       
C     SPECIAL ACTION WHEN N=1.                                          
C                                                                       
      IF (JU .NE. 0) GO TO 40                                           
         IF (CQY) QY(1) = Y(1)                                          
         IF (CQTY) QTY(1) = Y(1)                                        
         IF (CXB) XB(1) = Y(1)                                          
         IF (.NOT.CB) GO TO 30                                          
            IF (X(1,1) .NE. 0.0D0) GO TO 10                             
               INFO = 1                                                 
            GO TO 20                                                    
   10       CONTINUE                                                    
               B(1) = Y(1)/X(1,1)                                       
   20       CONTINUE                                                    
   30    CONTINUE                                                       
         IF (CR) RSD(1) = 0.0D0                                         
      GO TO 250                                                         
   40 CONTINUE                                                          
C                                                                       
C        SET UP TO COMPUTE QY OR QTY.                                   
C                                                                       
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)                                
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)                              
         IF (.NOT.CQY) GO TO 70                                         
C                                                                       
C           COMPUTE QY.                                                 
C                                                                       
            DO 60 JJ = 1, JU                                            
               J = JU - JJ + 1                                          
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50                        
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)              
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)                  
                  X(J,J) = TEMP                                         
   50          CONTINUE                                                 
   60       CONTINUE                                                    
   70    CONTINUE                                                       
         IF (.NOT.CQTY) GO TO 100                                       
C                                                                       
C           COMPUTE TRANS(Q)*Y.                                         
C                                                                       
            DO 90 J = 1, JU                                             
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80                        
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)             
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)                 
                  X(J,J) = TEMP                                         
   80          CONTINUE                                                 
   90       CONTINUE                                                    
  100    CONTINUE                                                       
C                                                                       
C        SET UP TO COMPUTE B, RSD, OR XB.                               
C                                                                       
         IF (CB) CALL DCOPY(K,QTY,1,B,1)                                
         KP1 = K + 1                                                    
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)                              
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)   
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120                        
            DO 110 I = KP1, N                                           
               XB(I) = 0.0D0                                            
  110       CONTINUE                                                    
  120    CONTINUE                                                       
         IF (.NOT.CR) GO TO 140                                         
            DO 130 I = 1, K                                             
               RSD(I) = 0.0D0                                           
  130       CONTINUE                                                    
  140    CONTINUE                                                       
         IF (.NOT.CB) GO TO 190                                         
C                                                                       
C           COMPUTE B.                                                  
C                                                                       
            DO 170 JJ = 1, K                                            
               J = K - JJ + 1                                           
               IF (X(J,J) .NE. 0.0D0) GO TO 150                         
                  INFO = J                                              
C           ......EXIT                                                  
                  GO TO 180                                             
  150          CONTINUE                                                 
               B(J) = B(J)/X(J,J)                                       
               IF (J .EQ. 1) GO TO 160                                  
                  T = -B(J)                                             
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)                        
  160          CONTINUE                                                 
  170       CONTINUE                                                    
  180       CONTINUE                                                    
  190    CONTINUE                                                       
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240                          
C                                                                       
C           COMPUTE RSD OR XB AS REQUIRED.                              
C                                                                       
            DO 230 JJ = 1, JU                                           
               J = JU - JJ + 1                                          
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220                       
                  TEMP = X(J,J)                                         
                  X(J,J) = QRAUX(J)                                     
                  IF (.NOT.CR) GO TO 200                                
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)          
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)              
  200             CONTINUE                                              
                  IF (.NOT.CXB) GO TO 210                               
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)           
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)               
  210             CONTINUE                                              
                  X(J,J) = TEMP                                         
  220          CONTINUE                                                 
  230       CONTINUE                                                    
  240    CONTINUE                                                       
  250 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE DRADB2(IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADB2
C***REFER TO  DRFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADB2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  DRADB2
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE DRADB3(IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADB3
C***REFER TO  DRFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADB3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  DRADB3
      TAUR = -.5D0
      TAUI = .5D0*SQRT(3.0D0)
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE DRADB4(IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADB4
C***REFER TO  DRFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADB4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  DRADB4
      SQRT2 = SQRT(2.0D0)
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE DRADB5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADB5
C***REFER TO  DRFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADB5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  DRADB5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(.1D0*PI)
      TI11 = SIN(.4D0*PI)
      TR12 = -SIN(.3D0*PI)
      TI12 = SIN(.2D0*PI)
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE DRADBG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADBG
C***REFER TO  DRFFTB
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADBG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
C***FIRST EXECUTABLE STATEMENT  DRADBG
      TPI = 8.0D0*ATAN(1.0D0)
      ARG = TPI/DBLE(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
CDIR$ IVDEP
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
CDIR$ IVDEP
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.0D0
      AI1 = 0.0D0
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
CDIR$ IVDEP
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
CDIR$ IVDEP
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
      SUBROUTINE DRADF2(IDO,L1,CC,CH,WA1)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADF2
C***REFER TO  DRFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADF2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           ,
     1                WA1(*)
C***FIRST EXECUTABLE STATEMENT  DRADF2
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE DRADF3(IDO,L1,CC,CH,WA1,WA2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADF3
C***REFER TO  DRFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADF3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
C***FIRST EXECUTABLE STATEMENT  DRADF3
      TAUR = -0.5D0
      TAUI = 0.5D0*SQRT(3.0D0)
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE DRADF4(IDO,L1,CC,CH,WA1,WA2,WA3)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADF4
C***REFER TO  DRFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADF4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
C***FIRST EXECUTABLE STATEMENT  DRADF4
      HSQT2 = 0.5D0*SQRT(2.0D0)
      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 111
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      GO TO 110
  111 DO 109 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 108 K=1,L1
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  108    CONTINUE
  109 CONTINUE
  110 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE DRADF5(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADF5
C***REFER TO  DRFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADF5
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
C***FIRST EXECUTABLE STATEMENT  DRADF5
      PI = 4.0D0*ATAN(1.0D0)
      TR11 = SIN(0.1D0*PI)
      TI11 = SIN(0.4D0*PI)
      TR12 = -SIN(0.3D0*PI)
      TI12 = SIN(0.2D0*PI)
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
      SUBROUTINE DRADFG(IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRADFG
C***REFER TO  DRFFTF
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRADFG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
C***FIRST EXECUTABLE STATEMENT  DRADFG
      TPI = 8.0D0*ATAN(1.0D0)
      ARG = TPI/DBLE(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
CDIR$ IVDEP
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
CDIR$ IVDEP
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.0D0
      AI1 = 0.0D0
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
CDIR$ IVDEP
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE DRFFTB(N,R,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFFTB
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  890413   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Backward transform of a dp coefficient array.
C***DESCRIPTION
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1989
C
C  Subroutine DRFFTB computes the dp perodic sequence from its
C  Fourier coefficients (Fourier synthesis).  The transform is defined
C  below at output parameter R.
C
C  Input Parameters
C
C  N       the length of the array R to be transformed.  The method
C          is most efficient when N is a product of small primes.
C          N may change so long as different work arrays are provided.
C
C  R       a dp array of length N which contains the sequence
C          to be transformed
C
C  WSAVE   a dp work array which must be dimensioned at least 2*N+15
C          in the program that calls DRFFTB.  The WSAVE array must be
C          initialized by calling subroutine DRFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          The same WSAVE array can be used by DRFFTF and DRFFTB.
C
C
C  Output Parameters
C
C  R       For N even and For I = 1,...,N
C
C               R(I) = R(1)+(-1)**(I-1)*R(N)
C
C                    plus the sum from K=2 to K=N/2 of
C
C                     2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                    -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C          For N odd and For I = 1,...,N
C
C               R(I) = R(1) plus the sum from K=2 to K=(N+1)/2 of
C
C                    2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                   -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C   *****  Note:
C               This transform is unnormalized since a call of DRFFTF
C               followed by a call of DRFFTB will multiply the input
C               sequence by N.
C
C  WSAVE   contains results which must not be destroyed between
C          calls of DRFFTB or DRFFTF.
C
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DRFTB1
C***END PROLOGUE  DRFFTB
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       R(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DRFFTB
      IF (N .EQ. 1) RETURN
      CALL DRFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE DRFFTF(N,R,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFFTF
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  890413   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Forward transform of a dp, periodic sequence.
C***DESCRIPTION
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C
C  Subroutine DRFFTF computes the Fourier coefficients of a dp
C  perodic sequence (Fourier analysis).  The transform is defined
C  below at output parameter R.
C
C  Input Parameters
C
C  N       the length of the array R to be transformed.  The method
C          is most efficient when N is a product of small primes.
C          N may change so long as different work arrays are provided
C
C  R       a dp array of length N which contains the sequence
C          to be transformed
C
C  WSAVE   a dp work array which must be dimensioned at least 2*N+15
C          in the program that calls DRFFTF.  The WSAVE array must be
C          initialized by calling subroutine DRFFTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C          the same WSAVE array can be used by DRFFTF and DRFFTB.
C
C
C  Output Parameters
C
C  R       R(1) = the sum from I=1 to I=N of R(I)
C
C          If N is even set L = N/2; if N is odd set L = (N+1)/2
C
C            then for K = 2,...,L
C
C               R(2*K-2) = the sum from I = 1 to I = N of
C
C                    R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C               R(2*K-1) = the sum from I = 1 to I = N of
C
C                   -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C          If N is even
C
C               R(N) = the sum from I = 1 to I = N of
C
C                    (-1)**(I-1)*R(I)
C
C   *****  Note:
C               This transform is unnormalized since a call of DRFFTF
C               followed by a call of DRFFTB will multiply the input
C               sequence by N.
C
C  WSAVE   contains results which must not be destroyed between
C          calls of DRFFTF or DRFFTB.
C
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DRFTF1
C***END PROLOGUE  DRFFTF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       R(*)       ,WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DRFFTF
      IF (N .EQ. 1) RETURN
      CALL DRFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE DRFFTI(N,WSAVE)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFFTI
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860115   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FOURIER TRANSFORM
C***AUTHOR  SWARZTRAUBER, P. N., (NCAR)
C***PURPOSE  Initialize for DRFFTF and DRFFTB.
C***DESCRIPTION
C
C           From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1989
C
C  Subroutine DRFFTI initializes the array WSAVE which is used in
C  both DRFFTF and DRFFTB.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.
C
C  Output Parameter
C
C  WSAVE   a DP work array which must be dimensioned at least 2*N+15.
C          The same work array can be used for both DRFFTF and DRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of DRFFTF or DRFFTB.
C
C  *********************************************************************
C  *                                                                   *
C  *   SUBPROGRAM REVISION HISTORY                                     *
C  *                                                                   *
C  *   06/01/79  -  Original version by Paul Swarztrauber.             *
C  *                Distributed by NCAR (ref. 1).                      *
C  *   04/01/83  -  SLATEC Common Math Library Subcommittee.           *
C  *                Modified to use SLATEC library source file format. *
C  *                Distributed in the SLATEC library (ref. 2).        *
C  *   01/15/86  -  Ron Boisvert, National Bureau of Standards.        *
C  *                Modified to convert to portable Fortran 77.        *
C  *                                                                   *
C  *   The changes introduced in the most recent modification are      *
C  *                                                                   *
C  *   (a) Dummy array size declarations (1) changed to (*)            *
C  *   (b) References to intrinsic function FLOAT changed to DP        *
C  *   (c) Mathematical constants previously coded in DATA state-      *
C  *       ments now computed at runtime using Fortran intrinsic       *
C  *       functions.  The affected variables are                      *
C  *                                                                   *
C  *          PI      SQRT2   SQRT3   TAUR    TR11    TR12             *
C  *          PIH     TSQRT2          TAUI    TI11    TI12             *
C  *          TPI     HSQT2                                            *
C  *                                                                   *
C  *   References                                                      *
C  *                                                                   *
C  *   1. P.N. Swarztrauber, Vectorizing the FFTs, in Parallel         *
C  *      Computations (G. Rodrigue, ed.), Academic Press, 1982,       *
C  *      pp. 51-83.                                                   *
C  *   2. B.L. Buzbee, The SLATEC Common Math Library, in Sources      *
C  *      and Development of Mathematical Software (W. Cowell, ed.),   *
C  *      Prentice-Hall, 1984, pp. 302-318.                            *
C  *                                                                   *
C  *********************************************************************
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DRFTI1
C***END PROLOGUE  DRFFTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  DRFFTI
      IF (N .EQ. 1) RETURN
      CALL DRFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      SUBROUTINE DRFTB1(N,C,CH,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFTB1
C***REFER TO  DRFFTB
C***ROUTINES CALLED  DRADB2,DRADB3,DRADB4,DRADB5,DRADBG
C***END PROLOGUE  DRFTB1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
C***FIRST EXECUTABLE STATEMENT  DRFTB1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL DRADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DRADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DRADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DRADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL DRADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DRADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL DRADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DRADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DRADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DRADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
      SUBROUTINE DRFTF1(N,C,CH,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFTF1
C***REFER TO  DRFTF
C***ROUTINES CALLED  DRADF2,DRADF3,DRADF4,DRADF5,DRADFG
C***END PROLOGUE  DRFTF1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
C***FIRST EXECUTABLE STATEMENT  DRFTF1
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL DRADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL DRADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL DRADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL DRADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL DRADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL DRADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL DRADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL DRADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL DRADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL DRADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
      SUBROUTINE DRFTI1(N,WA,IFAC)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRFTI1
C***REFER TO  DRFFTI
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRFTI1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  DRFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.D0*ATAN(1.D0)
      ARGH = TPI/N
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = LD*ARGH
            FI = 0.D0
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.D0
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DRNOR() 

c*********************************************************************72
c
C***BEGIN PROLOGUE  DRNOR
C***DATE WRITTEN   810915 (YYMMDD)
C***REVISION DATE  870419 (YYMMDD)
C***CATEGORY NO.  L6A14
C***KEYWORDS  RANDOM NUMBERS, NORMAL DEVIATES
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  GENERATES NORMAL RANDOM NUMBERS, WITH MEAN ZERO AND
C             UNIT STANDARD DEVIATION, OFTEN DENOTED N(0,1).
C
C***DESCRIPTION
C
C       DRNOR GENERATES NORMAL RANDOM NUMBERS WITH ZERO MEAN AND
C       UNIT STANDARD DEVIATION, OFTEN DENOTED N(0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C   USE
C       FIRST TIME....
C                   Z = DSTART(ISEED)
C                     HERE ISEED IS ANY  N O N - Z E R O  INTEGER.
C                     THIS CAUSES INITIALIZATION OF THE PROGRAM.
C                     DSTART RETURNS A DOUBLE PRECISION ECHO OF ISEED.
C
C       SUBSEQUENT TIMES...
C                   Z = DRNOR()
C                     CAUSES THE NEXT DOUBLE PRECISION RANDOM NUMBER
C                           TO BE RETURNED AS Z.
C
C.....................................................................
C                 TYPICAL USAGE
C
C                    DOUBLE PRECISION DSTART,DRNOR,Z
C                    INTEGER ISEED,I
C                    ISEED = 305
C                    Z = DSTART(ISEED)
C                    DO 1 I = 1,10
C                       Z = DRNOR()
C                       WRITE(*,'(1X,D20.15)') Z
C                 1  CONTINUE 
C                    END
C
C
C***REFERENCES  MARSAGLIA & TSANG, "A FAST, EASILY IMPLEMENTED
C                 METHOD FOR SAMPLING FROM DECREASING OR
C                 SYMMETRIC UNIMODAL DENSITY FUNCTIONS",
C                 PUBLISHED IN SIAM J SISC, JUNE 1984.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DRNOR
      DOUBLE PRECISION AA,B,C,C1,C2,PC,X,Y,XN,V(65),DSTART,U(17),
     *S,T,UN,VNI
      INTEGER J,IA,IB,IC,II,JJ,ID,III,JJJ,L
      SAVE U,II,JJ
C
      DATA AA,B,C/0.123758602991705622657860318807D+02
     +,0.487899177760378968003825536710D+00
     +,0.126770580788654778410032042685D+02/
      DATA C1,C2,PC,XN/0.9689279D0,1.301198D0
     +,0.195830333955546914251231662871D-01
     +,0.277699426966287548981739308903D+01/
      DATA (V(L), L=1,18)
     +/0.340945028703977890881422029377D+00
     +,0.457314591866933536170242641762D+00
     +,0.539779281611666939596931167481D+00
     +,0.606242679653048904661174385559D+00
     +,0.663169062764524854743428299352D+00
     +,0.713697459056025891183276150202D+00
     +,0.759612474933920604605610034675D+00
     +,0.802035600355531312751497342081D+00
     +,0.841722667978955421276418428136D+00
     +,0.879210223208313639290346470191D+00
     +,0.914894804386750587541168254518D+00
     +,0.949079113753090250631877133376D+00
     +,0.982000481239888200218207508382D+00
     +,0.101384923802994173461911276018D+01
     +,0.104478103674017351825847605485D+01
     +,0.107492538202855351339149779813D+01
     +,0.110439170226812581109973656162D+01
     +,0.113327377624394079212251428682D+01/
      DATA (V(L), L=19,37)
     +/0.116165303013393172842858957666D+01
     +,0.118960104083873798956793871425D+01
     +,0.121718147070087121538258873613D+01
     +,0.124445158789824683238161436879D+01
     +,0.127146348057211969375402099579D+01
     +,0.129826504188319751190626947962D+01
     +,0.132490078218086109537654808436D+01
     +,0.135141250993337129690631764473D+01
     +,0.137783991287001181384096757263D+01
     +,0.140422106355997540689484486002D+01
     +,0.143059286850269131403410180874D+01
     +,0.145699147613767157824869156623D+01
     +,0.148345265660321931119703498108D+01
     +,0.151001216431851991531882612256D+01
     +,0.153670609335952099134607533122D+01
     +,0.156357123503769104042967185962D+01
     +,0.159064544701425352365935513885D+01
     +,0.161796804367444698360816323707D+01
     +,0.164558021836908161542865488149D+01/
      DATA (V(L), L=38,55)
     +/0.167352550956703867146009214486D+01
     +,0.170185032506274055254533570699D+01
     +,0.173060454131778319060975251429D+01
     +,0.175984219903830120010946138955D+01
     +,0.178962232156657450014351836107D+01
     +,0.182000989013069176863519209140D+01
     +,0.185107702023027589942628767312D+01
     +,0.188290439759287281399927405628D+01
     +,0.191558305194303202395065401364D+01
     +,0.194921657491636060191700129156D+01
     +,0.198392392890568577258506733664D+01
     +,0.201984305290623555306662745946D+01
     +,0.205713555999009616804474181513D+01
     +,0.209599295624939161758467989658D+01
     +,0.213664502254438986524966622832D+01
     +,0.217937134039813565892460111431D+01
     +,0.222451750721601784110056845259D+01
     +,0.227251855485014779874266158018D+01/
      DATA (V(L), L=56,65)
     +/0.232393382009430256940425938218D+01
     +,0.237950077408282829688673722776D+01
     +,0.244022179797994340264326423618D+01
     +,0.250751170186531701106382130475D+01
     +,0.258346583522542956831304962942D+01
     +,0.267139159032083601869533973173D+01
     +,0.277699426966286466722522163764D+01
     +,0.277699426966286466722522163764D+01
     +,0.277699426966286466722522163764D+01
     +,0.277699426966286466722522163764D+01/
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH SEED 305.
      DATA U/
     *0.47196 09815 77884 75583 77897 24978D+00,
     *0.93032 34532 05669 57843 36396 32431D+00,
     *0.11016 17909 33730 83658 71279 44899D+00,
     *0.57150 19962 73139 51836 26387 57010D-01,
     *0.40246 75547 79738 26623 75385 03137D+00,
     *0.45118 19534 27459 48945 82794 56915D+00,
     *0.29607 61523 42721 10217 41299 54053D+00,
     *0.12820 21893 25888 11646 68796 22359D-01,
     *0.31427 46938 50973 60398 08532 59266D+00,
     *0.33552 13667 52294 93246 81635 94171D-02,
     *0.48868 50452 00439 37160 78503 67840D+00,
     *0.19547 04268 65656 75869 38606 13516D+00,
     *0.86416 27067 91773 55690 15993 26053D+00,
     *0.33550 59558 15259 20359 63811 70316D+00,
     *0.37719 02001 99058 08546 95264 70541D+00,
     *0.40078 03921 14818 31467 16765 25916D+00,
     *0.37422 42141 82207 46626 27503 07281D+00/
C
      DATA II,JJ / 17, 5 /
C
C***FIRST EXECUTABLE STATEMENT  DRNOR
C
C FAST PART...
C
C
C   BASIC GENERATOR IS FIBONACCI
C
      UN = U(II)-U(JJ)
      IF(UN.LT.0.0D0) UN = UN+1.
      U(II) = UN
C           U(II) AND UN ARE UNIFORM ON [0,1)
C           VNI IS UNIFORM ON [-1,1)
      VNI = UN + UN -1.
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C        INT(UN(II)*128) IN RANGE [0,127],  J IS IN RANGE [1,64]
      J = MOD(INT(U(II)*128),64)+1
C        PICK SIGN AS VNI IS POSITIVE OR NEGATIVE 
      DRNOR = VNI*V(J+1)
      IF(ABS(DRNOR).LE.V(J))RETURN
C
C SLOW PART; AA IS A*F(0)
      X = (ABS(DRNOR)-V(J))/(V(J+1)-V(J))
C          Y IS UNIFORM ON [0,1)
      Y = U(II)-U(JJ)
      IF(Y.LT.0.0D0) Y = Y+1. 
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C
      S = X+Y
      IF(S.GT.C2)GO TO 11
      IF(S.LE.C1)RETURN
      IF(Y.GT.C-AA*EXP(-.5D0*(B-B*X)**2))GO TO 11 
      IF(EXP(-.5D0*V(J+1)**2)+Y*PC/V(J+1).LE.EXP(-.5D0*DRNOR**2))
     *RETURN
C
C TAIL PART; .36010157... IS 1.0D0/XN
C       Y IS UNIFORM ON [0,1) 
   22 Y = U(II)-U(JJ)
      IF(Y.LE.0.0D0) Y = Y+1. 
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
C
      X = 0.360101571301190680192994239651D+00*LOG(Y)
C       Y IS UNIFORM ON [0,1) 
      Y = U(II)-U(JJ)
      IF(Y.LE.0.0D0) Y = Y+1. 
      U(II) = Y
      II = II-1
      IF(II.EQ.0)II = 17
      JJ = JJ-1
      IF(JJ.EQ.0)JJ = 17
      IF( -2.0D0*LOG(Y).LE.X**2 )GO TO 22
      DRNOR = SIGN(XN-X,DRNOR)
      RETURN
   11 DRNOR = SIGN(B-B*X,DRNOR)
      RETURN
C
C
C  FILL
      ENTRY DSTART(ISEED)
      IF(ISEED.NE.0) THEN
C
C          SET UP ...
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED
C
        II = 17
        JJ = 5
        IA = MOD(ABS(ISEED),32707)
        IB = 1111
        IC = 1947
        DO 2 III = 1,17
          S = 0.0D0 
          T = 0.5D0 
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER 95 BITS, ENOUGH FOR MOST MACHINES
C                   IN DOUBLE PRECISION.
          DO 3 JJJ = 1,95
                  ID = IC-IA
                  IF(ID.GE.0)GOTO 4
                  ID = ID+32707
                  S = S+T
    4             IA = IB
                  IB = IC
                  IC = ID
    3     T = 0.5D0*T
    2     U(III) = S
      ENDIF
C       RETURN FLOATING ECHO OF ISEED
      DSTART=ISEED
      RETURN
      END
      subroutine drot ( n, dx, incx, dy, incy, c, s )

c*********************************************************************72
c
cc DROT applies a plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, double precision X(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Input, double precision C, S, parameters (presumably the cosine and
c    sine of some angle) that define a plane rotation.
c
      implicit none

      double precision dx(*),dy(*),dtemp,c,s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
      end do
      return
      end
      subroutine drotg ( da, db, c, s )

c*********************************************************************72
c
cc DROTG constructs a Givens plane rotation.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    Given values A and B, this routine computes
c
c    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
c          = sign ( B ) if abs ( A ) <= abs ( B );
c
c    R     = SIGMA * ( A * A + B * B );
c
c    C = A / R if R is not 0
c      = 1     if R is 0;
c
c    S = B / R if R is not 0,
c        0     if R is 0.
c
c    The computed numbers then satisfy the equation
c
c    (  C  S ) ( A ) = ( R )
c    ( -S  C ) ( B ) = ( 0 )
c
c    The routine also computes
c
c    Z = S     if abs ( A ) > abs ( B ),
c      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
c      = 1     if C is 0.
c
c    The single value Z encodes C and S, and hence the rotation:
c
c    If Z = 1, set C = 0 and S = 1;
c    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
c    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input/output, double precision SA, SB.  On input, SA and SB are the values
c    A and B.  On output, SA is overwritten with R, and SB is
c    overwritten with Z.
c
c    Output, double precision C, S, the cosine and sine of the
c    Givens rotation.
c
      implicit none

      double precision da,db,c,s,roe,scale,r,z

      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)

      if( scale .eq. 0.0d0 ) then
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
      else
        r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
        r = dsign(1.0d0,roe)*r
        c = da/r
        s = db/r
        z = 1.0d0
        if( dabs(da) .gt. dabs(db) ) z = s
        if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
      end if

      da = r
      db = z

      return
      end
      SUBROUTINE DSCAL(N,DA,DX,INCX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SSCAL-S DSCAL-D CSCAL-C),
C             LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSCAL
C
      DOUBLE PRECISION DA,DX(1)
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DSVDC
C***DATE WRITTEN   790319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D6
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,
C             SINGULAR VALUE DECOMPOSITION
C***AUTHOR  STEWART, G. W., (U. OF MARYLAND)
C***PURPOSE  Perform the singular value decomposition of a d.p. NXP
C            matrix.
C***DESCRIPTION
C
C     DSVDC is a subroutine to reduce a double precision NxP matrix X
C     by orthogonal transformations U and V to diagonal form.  The
C     diagonal elements S(I) are the singular values of X.  The
C     columns of U are the corresponding left singular vectors,
C     and the columns of V the right singular vectors.
C
C     On Entry
C
C         X         DOUBLE PRECISION(LDX,P), where LDX .GE. N.
C                   X contains the matrix whose singular value
C                   decomposition is to be computed.  X is
C                   destroyed by DSVDC.
C
C         LDX       INTEGER.
C                   LDX is the leading dimension of the array X.
C
C         N         INTEGER.
C                   N is the number of columns of the matrix X.
C
C         P         INTEGER.
C                   P is the number of rows of the matrix X.
C
C         LDU       INTEGER.
C                   LDU is the leading dimension of the array U.
C                   (See below).
C
C         LDV       INTEGER.
C                   LDV is the leading dimension of the array V.
C                   (See below).
C
C         WORK      DOUBLE PRECISION(N).
C                   WORK is a scratch array.
C
C         JOB       INTEGER.
C                   JOB controls the computation of the singular
C                   vectors.  It has the decimal expansion AB
C                   with the following meaning
C
C                        A .EQ. 0    do not compute the left singular
C                                  vectors.
C                        A .EQ. 1    return the N left singular vectors
C                                  in U.
C                        A .GE. 2    return the first MIN(N,P) singular
C                                  vectors in U.
C                        B .EQ. 0    do not compute the right singular
C                                  vectors.
C                        B .EQ. 1    return the right singular vectors
C                                  in V.
C
C     On Return
C
C         S         DOUBLE PRECISION(MM), where MM=MIN(N+1,P).
C                   The first MIN(N,P) entries of S contain the
C                   singular values of X arranged in descending
C                   order of magnitude.
C
C         E         DOUBLE PRECISION(P).
C                   E ordinarily contains zeros.  However see the
C                   discussion of INFO for exceptions.
C
C         U         DOUBLE PRECISION(LDU,K), where LDU .GE. N.
C                   If JOBA .EQ. 1, then K .EQ. N.
C                   If JOBA .GE. 2, then K .EQ. MIN(N,P).
C                   U contains the matrix of right singular vectors.
C                   U is not referenced if JOBA .EQ. 0.  If N .LE. P
C                   or if JOBA .EQ. 2, then U may be identified with X
C                   in the subroutine call.
C
C         V         DOUBLE PRECISION(LDV,P), where LDV .GE. P.
C                   V contains the matrix of right singular vectors.
C                   V is not referenced if JOB .EQ. 0.  If P .LE. N,
C                   then V may be identified with X in the
C                   subroutine call.
C
C         INFO      INTEGER.
C                   The singular values (and their corresponding
C                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
C                   are correct (here M=MIN(N,P)).  Thus if
C                   INFO .EQ. 0, all the singular values and their
C                   vectors are correct.  In any event, the matrix
C                   B = TRANS(U)*X*V is the bidiagonal matrix
C                   with the elements of S on its diagonal and the
C                   elements of E on its super-diagonal (TRANS(U)
C                   is the transpose of U).  Thus the singular
C                   values of X and B are the same.
C
C     LINPACK.  This version dated 03/19/79 .
C     G. W. Stewart, University of Maryland, Argonne National Lab.
C
C     DSVDC uses the following functions and subprograms.
C
C     External DROT
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
C     Fortran DABS,DMAX1,MAX0,MIN0,MOD,DSQRT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT,DNRM2,DROT,DROTG,DSCAL,DSWAP
C***END PROLOGUE  DSVDC
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO
      DOUBLE PRECISION X(LDX,1),S(1),E(1),U(LDU,1),V(LDV,1),WORK(1)
C
C
      INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT,
     1        MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION B,C,CS,EL,EMM1,F,G,DNRM2,SCALE,SHIFT,SL,SM,SN,
     1                 SMM1,T1,TEST,ZTEST
      LOGICAL WANTU,WANTV
C
C     SET THE MAXIMUM NUMBER OF ITERATIONS.
C
C***FIRST EXECUTABLE STATEMENT  DSVDC
      MAXIT = 30
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU .GT. 1) NCU = MIN0(N,P)
      IF (JOBU .NE. 0) WANTU = .TRUE.
      IF (MOD(JOB,10) .NE. 0) WANTV = .TRUE.
C
C     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
C     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
C
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU .LT. 1) GO TO 170
      DO 160 L = 1, LU
         LP1 = L + 1
         IF (L .GT. NCT) GO TO 20
C
C           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
C           PLACE THE L-TH DIAGONAL IN S(L).
C
            S(L) = DNRM2(N-L+1,X(L,L),1)
            IF (S(L) .EQ. 0.0D0) GO TO 10
               IF (X(L,L) .NE. 0.0D0) S(L) = DSIGN(S(L),X(L,L))
               CALL DSCAL(N-L+1,1.0D0/S(L),X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
   10       CONTINUE
            S(L) = -S(L)
   20    CONTINUE
         IF (P .LT. LP1) GO TO 50
         DO 40 J = LP1, P
            IF (L .GT. NCT) GO TO 30
            IF (S(L) .EQ. 0.0D0) GO TO 30
C
C              APPLY THE TRANSFORMATION.
C
               T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
               CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
C
C           PLACE THE L-TH ROW OF X INTO  E FOR THE
C           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
C
            E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
         IF (.NOT.WANTU .OR. L .GT. NCT) GO TO 70
C
C           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
C           MULTIPLICATION.
C
            DO 60 I = L, N
               U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
         IF (L .GT. NRT) GO TO 150
C
C           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
C           L-TH SUPER-DIAGONAL IN E(L).
C
            E(L) = DNRM2(P-L,E(LP1),1)
            IF (E(L) .EQ. 0.0D0) GO TO 80
               IF (E(LP1) .NE. 0.0D0) E(L) = DSIGN(E(L),E(LP1))
               CALL DSCAL(P-L,1.0D0/E(L),E(LP1),1)
               E(LP1) = 1.0D0 + E(LP1)
   80       CONTINUE
            E(L) = -E(L)
            IF (LP1 .GT. N .OR. E(L) .EQ. 0.0D0) GO TO 120
C
C              APPLY THE TRANSFORMATION.
C
               DO 90 I = LP1, N
                  WORK(I) = 0.0D0
   90          CONTINUE
               DO 100 J = LP1, P
                  CALL DAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
               DO 110 J = LP1, P
                  CALL DAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
            IF (.NOT.WANTV) GO TO 140
C
C              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
C              BACK MULTIPLICATION.
C
               DO 130 I = LP1, P
                  V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
C
C     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
C
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT .LT. P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N .LT. M) S(M) = 0.0D0
      IF (NRTP1 .LT. M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0D0
C
C     IF REQUIRED, GENERATE U.
C
      IF (.NOT.WANTU) GO TO 300
         IF (NCU .LT. NCTP1) GO TO 200
         DO 190 J = NCTP1, NCU
            DO 180 I = 1, N
               U(I,J) = 0.0D0
  180       CONTINUE
            U(J,J) = 1.0D0
  190    CONTINUE
  200    CONTINUE
         IF (NCT .LT. 1) GO TO 290
         DO 280 LL = 1, NCT
            L = NCT - LL + 1
            IF (S(L) .EQ. 0.0D0) GO TO 250
               LP1 = L + 1
               IF (NCU .LT. LP1) GO TO 220
               DO 210 J = LP1, NCU
                  T = -DDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
                  CALL DAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
               CALL DSCAL(N-L+1,-1.0D0,U(L,L),1)
               U(L,L) = 1.0D0 + U(L,L)
               LM1 = L - 1
               IF (LM1 .LT. 1) GO TO 240
               DO 230 I = 1, LM1
                  U(I,L) = 0.0D0
  230          CONTINUE
  240          CONTINUE
            GO TO 270
  250       CONTINUE
               DO 260 I = 1, N
                  U(I,L) = 0.0D0
  260          CONTINUE
               U(L,L) = 1.0D0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C     IF IT IS REQUIRED, GENERATE V.
C
      IF (.NOT.WANTV) GO TO 350
         DO 340 LL = 1, P
            L = P - LL + 1
            LP1 = L + 1
            IF (L .GT. NRT) GO TO 320
            IF (E(L) .EQ. 0.0D0) GO TO 320
               DO 310 J = LP1, P
                  T = -DDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
                  CALL DAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
            DO 330 I = 1, P
               V(I,L) = 0.0D0
  330       CONTINUE
            V(L,L) = 1.0D0
  340    CONTINUE
  350 CONTINUE
C
C     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
C
      MM = M
      ITER = 0
  360 CONTINUE
C
C        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
C
C     ...EXIT
         IF (M .EQ. 0) GO TO 620
C
C        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
C        FLAG AND RETURN.
C
         IF (ITER .LT. MAXIT) GO TO 370
            INFO = M
C     ......EXIT
            GO TO 620
  370    CONTINUE
C
C        THIS SECTION OF THE PROGRAM INSPECTS FOR
C        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
C        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
C
C           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
C           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
C           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
C                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
C           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
C
         DO 390 LL = 1, M
            L = M - LL
C        ...EXIT
            IF (L .EQ. 0) GO TO 400
            TEST = DABS(S(L)) + DABS(S(L+1))
            ZTEST = TEST + DABS(E(L))
            IF (ZTEST .NE. TEST) GO TO 380
               E(L) = 0.0D0
C        ......EXIT
               GO TO 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
         IF (L .NE. M - 1) GO TO 410
            KASE = 4
         GO TO 480
  410    CONTINUE
            LP1 = L + 1
            MP1 = M + 1
            DO 430 LLS = LP1, MP1
               LS = M - LLS + LP1
C           ...EXIT
               IF (LS .EQ. L) GO TO 440
               TEST = 0.0D0
               IF (LS .NE. M) TEST = TEST + DABS(E(LS))
               IF (LS .NE. L + 1) TEST = TEST + DABS(E(LS-1))
               ZTEST = TEST + DABS(S(LS))
               IF (ZTEST .NE. TEST) GO TO 420
                  S(LS) = 0.0D0
C           ......EXIT
                  GO TO 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
            IF (LS .NE. L) GO TO 450
               KASE = 3
            GO TO 470
  450       CONTINUE
            IF (LS .NE. M) GO TO 460
               KASE = 1
            GO TO 470
  460       CONTINUE
               KASE = 2
               L = LS
  470       CONTINUE
  480    CONTINUE
         L = L + 1
C
C        PERFORM THE TASK INDICATED BY KASE.
C
         GO TO (490,520,540,570), KASE
C
C        DEFLATE NEGLIGIBLE S(M).
C
  490    CONTINUE
            MM1 = M - 1
            F = E(M-1)
            E(M-1) = 0.0D0
            DO 510 KK = L, MM1
               K = MM1 - KK + L
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               IF (K .EQ. L) GO TO 500
                  F = -SN*E(K-1)
                  E(K-1) = CS*E(K-1)
  500          CONTINUE
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
         GO TO 610
C
C        SPLIT AT NEGLIGIBLE S(L).
C
  520    CONTINUE
            F = E(L-1)
            E(L-1) = 0.0D0
            DO 530 K = L, M
               T1 = S(K)
               CALL DROTG(T1,F,CS,SN)
               S(K) = T1
               F = -SN*E(K)
               E(K) = CS*E(K)
               IF (WANTU) CALL DROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
         GO TO 610
C
C        PERFORM ONE QR STEP.
C
  540    CONTINUE
C
C           CALCULATE THE SHIFT.
C
            SCALE = DMAX1(DABS(S(M)),DABS(S(M-1)),DABS(E(M-1)),
     1                    DABS(S(L)),DABS(E(L)))
            SM = S(M)/SCALE
            SMM1 = S(M-1)/SCALE
            EMM1 = E(M-1)/SCALE
            SL = S(L)/SCALE
            EL = E(L)/SCALE
            B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0D0
            C = (SM*EMM1)**2
            SHIFT = 0.0D0
            IF (B .EQ. 0.0D0 .AND. C .EQ. 0.0D0) GO TO 550
               SHIFT = DSQRT(B**2+C)
               IF (B .LT. 0.0D0) SHIFT = -SHIFT
               SHIFT = C/(B + SHIFT)
  550       CONTINUE
            F = (SL + SM)*(SL - SM) - SHIFT
            G = SL*EL
C
C           CHASE ZEROS.
C
            MM1 = M - 1
            DO 560 K = L, MM1
               CALL DROTG(F,G,CS,SN)
               IF (K .NE. L) E(K-1) = F
               F = CS*S(K) + SN*E(K)
               E(K) = CS*E(K) - SN*S(K)
               G = SN*S(K+1)
               S(K+1) = CS*S(K+1)
               IF (WANTV) CALL DROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
               CALL DROTG(F,G,CS,SN)
               S(K) = F
               F = CS*E(K) + SN*S(K+1)
               S(K+1) = -SN*E(K) + CS*S(K+1)
               G = SN*E(K+1)
               E(K+1) = CS*E(K+1)
               IF (WANTU .AND. K .LT. N)
     1            CALL DROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
            E(M-1) = F
            ITER = ITER + 1
         GO TO 610
C
C        CONVERGENCE.
C
  570    CONTINUE
C
C           MAKE THE SINGULAR VALUE  POSITIVE.
C
            IF (S(L) .GE. 0.0D0) GO TO 580
               S(L) = -S(L)
               IF (WANTV) CALL DSCAL(P,-1.0D0,V(1,L),1)
  580       CONTINUE
C
C           ORDER THE SINGULAR VALUE.
C
  590       IF (L .EQ. MM) GO TO 600
C           ...EXIT
               IF (S(L) .GE. S(L+1)) GO TO 600
               T = S(L)
               S(L) = S(L+1)
               S(L+1) = T
               IF (WANTV .AND. L .LT. P)
     1            CALL DSWAP(P,V(1,L),1,V(1,L+1),1)
               IF (WANTU .AND. L .LT. N)
     1            CALL DSWAP(N,U(1,L),1,U(1,L+1),1)
               L = L + 1
            GO TO 590
  600       CONTINUE
            ITER = 0
            M = M - 1
  610    CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END
      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)

c*********************************************************************72
c
C***BEGIN PROLOGUE  DSWAP
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A5
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(SSWAP-S DSWAP-D CSWAP-C ISWAP-I),
C             INTERCHANGE,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Interchange d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSWAP
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP1,DTEMP2,DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DUNI()

c*********************************************************************72
c
C***BEGIN PROLOGUE  DUNI
C***DATE WRITTEN   880714 (YYMMDD)
C***REVISION DATE  880714 (YYMMDD)
C***CATEGORY NO.  L6A21
C***KEYWORDS  RANDOM NUMBERS, UNIFORM RANDOM NUMBERS
C***AUTHOR    KAHANER, DAVID, SCIENTIFIC COMPUTING DIVISION, NBS
C             MARSAGLIA, GEORGE, SUPERCOMPUTER RES. INST., FLORIDA ST. U.
C
C***PURPOSE  THIS ROUTINE GENERATES DOUBLE PRECISION UNIFORM
C             RANDOM NUMBERS ON [0,1)
C***DESCRIPTION
C        COMPUTES DOUBLE PRECISION UNIFORM NUMBERS ON [0,1).
C           FROM THE BOOK, "NUMERICAL METHODS AND SOFTWARE" BY
C                D. KAHANER, C. MOLER, S. NASH
C                PRENTICE HALL, 1988
C
C       USAGE: 
C              TO INITIALIZE THE GENERATOR
C                   USEED = DUSTAR(ISEED)
C               WHERE: ISEED IS ANY NONZERO INTEGER
C                  WILL RETURN FLOATING POINT VALUE OF ISEED.
C
C               SUBSEQUENTLY
C                       U = DUNI()
C                  WILL RETURN A REAL UNIFORM ON [0,1)
C
C                ONE INITIALIZATION IS NECESSARY, BUT ANY NUMBER OF EVALUATIONS 
C                  OF DUNI IN ANY ORDER, ARE ALLOWED.
C
C           NOTE: DEPENDING UPON THE VALUE OF K (SEE BELOW), THE OUTPUT
C                       OF DUNI MAY DIFFER FROM ONE MACHINE TO ANOTHER.
C
C           TYPICAL USAGE: 
C
C               DOUBLE PRECISION U,DUNI,DUSTAR,USEED
C               INTEGER ISEED 
CC                 SET SEED
C               ISEED = 305
C               USEED = DUSTAR(ISEED)
C               DO 1 I = 1,1000
C                   U = DUNI()
C             1 CONTINUE
CC                 NOTE: IF K=47 (THE DEFAULT, SEE BELOW) THE OUTPUT VALUE OF
CC                           U WILL BE 0.812053811384E-01...
C               WRITE(*,*) U
C               END 
C
C          NOTE ON PORTABILITY: USERS CAN CHOOSE TO RUN DUNI IN ITS DEFAULT
C               MODE (REQUIRING NO USER ACTION) WHICH WILL GENERATE THE SAME
C               SEQUENCE OF NUMBERS ON ANY COMPUTER SUPPORTING FLOATING POINT
C               NUMBERS WITH AT LEAST 47 BIT MANTISSAS, OR IN A MODE THAT
C               WILL GENERATE NUMBERS WITH A LONGER PERIOD ON COMPUTERS WITH
C               LARGER MANTISSAS.
C          TO EXERCISE THIS OPTION:  B E F O R E  INVOKING DUSTAR INSERT
C               THE INSTRUCTION        UBITS = DUNIB(K)      K >= 47
C               WHERE K IS THE NUMBER OF BITS IN THE MANTISSA OF YOUR FLOATING
C               POINT WORD (K=96 FOR CRAY, CYBER 205). DUNIB RETURNS THE
C               FLOATING POINT VALUE OF K THAT IT ACTUALLY USED.
C                    K INPUT AS .LE. 47, THEN UBITS=47.
C                    K INPUT AS .GT. 47, THEN UBITS=FLOAT(K)
C               IF K>47 THE SEQUENCE OF NUMBERS GENERATED BY DUNI MAY DIFFER
C               FROM ONE COMPUTER TO ANOTHER.
C
C
C
C***REFERENCES  MARSAGLIA G., "COMMENTS ON THE PERFECT UNIFORM RANDOM 
C                 NUMBER GENERATOR", UNPUBLISHED NOTES, WASH S. U.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE DUNI
      DOUBLE PRECISION CSAVE,CD,CM
      PARAMETER(
     *    CSAVE=0.9162596898123D+13/0.140737488355328D+15,
     *    CD=0.76543212345678D+14/0.140737488355328D+15,
     *    CM=0.140737488355213D+15/0.140737488355328D+15)
C                            2**47=0.140737488355328D+15
      DOUBLE PRECISION U(17),S,T,DUSTAR,C,DUNIB
      INTEGER I,J,II,JJ,K,KK,I1,J1,K1,L1,M1,ISEED 
C
      SAVE U,I,J,K,C
C      LOAD DATA ARRAY IN CASE USER FORGETS TO INITIALIZE.
C      THIS ARRAY IS THE RESULT OF CALLING DUNI 100000 TIMES
C         WITH ISEED=305 AND K=96.
      DATA U/
     *0.471960981577884755837789724978D+00,
     *0.930323453205669578433639632431D+00,
     *0.110161790933730836587127944899D+00,
     *0.571501996273139518362638757010D-01,
     *0.402467554779738266237538503137D+00,
     *0.451181953427459489458279456915D+00,
     *0.296076152342721102174129954053D+00,
     *0.128202189325888116466879622359D-01,
     *0.314274693850973603980853259266D+00,
     *0.335521366752294932468163594171D-02,
     *0.488685045200439371607850367840D+00,
     *0.195470426865656758693860613516D+00,
     *0.864162706791773556901599326053D+00,
     *0.335505955815259203596381170316D+00,
     *0.377190200199058085469526470541D+00,
     *0.400780392114818314671676525916D+00,
     *0.374224214182207466262750307281D+00/
      DATA I,J,K,C/17,5,47,CSAVE/
C
C   BASIC GENERATOR IS FIBONACCI
C
      DUNI = U(I)-U(J)
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      U(I) = DUNI
      I = I-1
      IF(I.EQ.0)I = 17
      J = J-1
      IF(J.EQ.0)J = 17
C
C   SECOND GENERATOR IS CONGRUENTIAL
C
      C = C-CD
      IF(C.LT.0.0D0) C=C+CM
C
C   COMBINATION GENERATOR
C
      DUNI = DUNI-C 
      IF(DUNI.LT.0.0D0)DUNI = DUNI+1.0D0
      RETURN
C
      ENTRY DUSTAR(ISEED)
C
C          SET UP ...
C          CONVERT ISEED TO FOUR SMALLISH POSITIVE INTEGERS.
C
        I1 = MOD(ABS(ISEED),177)+1
        J1 = MOD(ABS(ISEED),167)+1
        K1 = MOD(ABS(ISEED),157)+1
        L1 = MOD(ABS(ISEED),147)+1
C
C              GENERATE RANDOM BIT PATTERN IN ARRAY BASED ON GIVEN SEED.
C
        DO 2 II = 1,17
          S = 0.0D0 
          T = 0.5D0 
C             DO FOR EACH OF THE BITS OF MANTISSA OF WORD
C             LOOP  OVER K BITS, WHERE K IS DEFAULTED TO 47 BUT CAN
C               BE CHANGED BY USER CALL TO DUNIB(K)
          DO 3 JJ = 1,K
                  M1 = MOD(MOD(I1*J1,179)*K1,179) 
                  I1 = J1
                  J1 = K1
                  K1 = M1
                  L1 = MOD(53*L1+1,169) 
                  IF(MOD(L1*M1,64).GE.32)S=S+T
    3             T = 0.5D0*T 
    2   U(II) = S
        DUSTAR = FLOAT(ISEED) 
        RETURN
C
      ENTRY DUNIB(KK)
        IF(KK.LE.47)THEN
             K=47
        ELSE
             K=KK
        ENDIF
        DUNIB=FLOAT(K)
      END
      SUBROUTINE FDUMP

c*********************************************************************72
c
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      SUBROUTINE FORSLD(NR,N,A,X,B)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C SOLVE  AX=B  WHERE A IS LOWER TRIANGULAR MATRIX
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED)
C X(N)        <--  SOLUTION VECTOR
C B(N)         --> RIGHT-HAND SIDE VECTOR
C
C NOTE
C ----
C IF B IS NO LONGER REQUIRED BY CALLING ROUTINE,
C THEN VECTORS B AND X MAY SHARE THE SAME STORAGE.
C
      DIMENSION A(NR,1),X(N),B(N)
C
C SOLVE LX=B. (FOREWARD SOLVE)
C
      X(1)=B(1)/A(1,1)
      IF(N.EQ.1) RETURN
      DO 20 I=2,N
        SUM=0.0D0
        IM1=I-1
        DO 10 J=1,IM1
          SUM=SUM+A(I,J)*X(J)
   10   CONTINUE
        X(I)=(B(I)-SUM)/A(I,I)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE FSTCDD (N, X, FCN, SX, RNOISE, G)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PURPOSE
C -------
C FIND CENTRAL DIFFERENCE APPROXIMATION G TO THE FIRST DERIVATIVE
C (GRADIENT) OF THE FUNCTION DEFINED BY FCN AT THE POINT X.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X            --> POINT AT WHICH GRADIENT IS TO BE APPROXIMATED.
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION.
C SX           --> DIAGONAL SCALING MATRIX FOR X.
C RNOISE       --> RELATIVE NOISE IN FCN [F(X)].
C G           <--  CENTRAL DIFFERENCE APPROXIMATION TO GRADIENT.
C
C
      DIMENSION X(N)
      DIMENSION SX(N)
      DIMENSION G(N)
      EXTERNAL FCN
C
C FIND I TH  STEPSIZE, EVALUATE TWO NEIGHBORS IN DIRECTION OF I TH
C UNIT VECTOR, AND EVALUATE I TH  COMPONENT OF GRADIENT.
C
      THIRD = 1.0D0/3.0D0
      DO 10 I = 1, N
         STEPI = RNOISE**THIRD * MAX(ABS(X(I)), 1.0D0/SX(I))
         XTEMPI = X(I)
         X(I) = XTEMPI + STEPI
         CALL FCN (N, X, FPLUS)
         X(I) = XTEMPI - STEPI
         CALL FCN (N, X, FMINUS)
         X(I) = XTEMPI
         G(I) = (FPLUS - FMINUS)/(2.0D0*STEPI)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FSTFDD(NR,M,N,XPLS,FCN,FPLS,A,SX,RNOISE,FHAT,ICASE)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PURPOSE
C -------
C FIND FIRST ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" TO THE
C FIRST DERIVATIVE OF THE FUNCTION DEFINED BY THE SUBPROGRAM "FNAME"
C EVALUATED AT THE NEW ITERATE "XPLS".
C
C
C FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE:
C 1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION FUNCTION "FCN
C    ANALYTIC USER ROUTINE HAS BEEN SUPPLIED;
C 2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION
C    IF NO ANALYTIC USER ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT
C    ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND IF THE
C    OPTIMIZATION FUNCTION IS INEXPENSIVE TO EVALUATE
C
C NOTE
C ----
C _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE FUNCTION
C      (FCN).   FCN(X) # F: R(N)-->R(1)
C _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE FUNCTION
C      FCN(X) # F: R(N)-->R(N).
C _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO
C      FUNCTION, WHERE THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN"
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C M            --> NUMBER OF ROWS IN A
C N            --> NUMBER OF COLUMNS IN A; DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE:  X[K]
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C FPLS(M)      --> _M=1 (OPTIMIZATION) FUNCTION VALUE AT NEW ITERATE:
C                       FCN(XPLS)
C                  _M=N (OPTIMIZATION) VALUE OF FIRST DERIVATIVE
C                       (GRADIENT) GIVEN BY USER FUNCTION FCN
C                  _M=N (SYSTEMS)  FUNCTION VALUE OF ASSOCIATED
C                       MINIMIZATION FUNCTION
C A(NR,N)     <--  FINITE DIFFERENCE APPROXIMATION (SEE NOTE).  ONLY
C                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE RETURNED
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNOISE       --> RELATIVE NOISE IN FCN [F(X)]
C FHAT(M)      --> WORKSPACE
C ICASE        --> =1 OPTIMIZATION (GRADIENT)
C                  =2 SYSTEMS
C                  =3 OPTIMIZATION (HESSIAN)
C
C INTERNAL VARIABLES
C ------------------
C STEPSZ - STEPSIZE IN THE J-TH VARIABLE DIRECTION
C
      DIMENSION XPLS(N),FPLS(M)
      DIMENSION FHAT(M)
      DIMENSION SX(N)
      DIMENSION A(NR,1)
C
C FIND J-TH COLUMN OF A
C EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J)
C
      DO 30 J=1,N
        STEPSZ=SQRT(RNOISE)*MAX(ABS(XPLS(J)),1.D0/SX(J))
        XTMPJ=XPLS(J)
        XPLS(J)=XTMPJ+STEPSZ
        CALL FCN(N,XPLS,FHAT)
        XPLS(J)=XTMPJ
        DO 20 I=1,M
          A(I,J)=(FHAT(I)-FPLS(I))/STEPSZ
   20   CONTINUE
   30 CONTINUE
      IF(ICASE.NE.3) RETURN
C
C IF COMPUTING HESSIAN, A MUST BE SYMMETRIC
C
      IF(N.EQ.1) RETURN
      NM1=N-1
      DO 50 J=1,NM1
        JP1=J+1
        DO 40 I=JP1,M
          A(I,J)=(A(I,J)+A(J,I))/2.0D0
   40   CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE GRCHKD(N,X,FCN,F,G,TYPSIZ,SX,FSCALE,RNF,
     +     ANALTL,WRK1,MSG,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C CHECK ANALYTIC GRADIENT AGAINST ESTIMATED GRADIENT
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ESTIMATE TO A ROOT OF FCN
C FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C                       FCN:  R(N) --> R(1)
C F            --> FUNCTION VALUE:  FCN(X)
C G(N)         --> GRADIENT:  G(X)
C TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X
C SX(N)        --> DIAGONAL SCALING MATRIX:  SX(I)=1./TYPSIZ(I)
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN
C RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN
C ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND
C                  ANALYTICAL GRADIENTS
C WRK1(N)      --> WORKSPACE
C MSG         <--  MESSAGE OR ERROR CODE
C                    ON OUTPUT: =-21, PROBABLE CODING ERROR OF GRADIENT
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),G(N)
      DIMENSION SX(N),TYPSIZ(N)
      DIMENSION WRK1(N)
      EXTERNAL FCN
C
C COMPUTE FIRST ORDER FINITE DIFFERENCE GRADIENT AND COMPARE TO
C ANALYTIC GRADIENT.
C
      CALL FSTFDD(1,1,N,X,FCN,F,WRK1,SX,RNF,WRK,1)
      KER=0
      DO 5 I=1,N
        GS=MAX(ABS(F),FSCALE)/MAX(ABS(X(I)),TYPSIZ(I))
        IF(ABS(G(I)-WRK1(I)).GT.MAX(ABS(G(I)),GS)*ANALTL) KER=1
    5 CONTINUE
      IF(KER.EQ.0) GO TO 20
        WRITE(IPR,901)
        WRITE(IPR,902) (I,G(I),WRK1(I),I=1,N)
        MSG=-21
   20 CONTINUE
      RETURN
  901 FORMAT(47H0GRCHKD    PROBABLE ERROR IN CODING OF ANALYTIC,
     +       19H GRADIENT FUNCTION./
     +       16H GRCHKD     COMP,12X,8HANALYTIC,12X,8HESTIMATE)
  902 FORMAT(11H GRCHKD    ,I5,3X,E20.13,3X,E20.13)
      END
      SUBROUTINE HOKDRD(NR,N,X,F,G,A,UDIAG,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,AMU,DLTP,PHI,PHIP0,
     +     SC,XPLSP,WRK0,EPSM,ITNCNT,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE (XPLS) BY THE MORE-HEBDON METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN LOWER
C                  TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN IN UPPER TRIANGULAR PART AND UDIAG.
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C P(N)         --> NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C IRETCD      <--  RETURN CODE
C                    =0 SATISFACTORY XPLS FOUND
C                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY
C                       DISTINCT FROM X
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C DLTP        <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C SC(N)        --> WORKSPACE
C XPLSP(N)     --> WORKSPACE
C WRK0(N)      --> WORKSPACE
C EPSM         --> MACHINE EPSILON
C ITNCNT       --> ITERATION COUNT
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),G(N),P(N),XPLS(N),SX(N)
      DIMENSION A(NR,1),UDIAG(N)
      DIMENSION SC(N),XPLSP(N),WRK0(N)
      LOGICAL MXTAKE,NWTAKE
      LOGICAL FSTIME
      EXTERNAL FCN
C
      IRETCD=4
      FSTIME=.TRUE.
      TMP=0.D0
      DO 5 I=1,N
        TMP=TMP+SX(I)*SX(I)*P(I)*P(I)
    5 CONTINUE
      RNWTLN=SQRT(TMP)
C$    WRITE(IPR,954) RNWTLN
C
      IF(ITNCNT.GT.1) GO TO 100
C     IF(ITNCNT.EQ.1)
C     THEN
        AMU=0.D0
C
C       IF FIRST ITERATION AND TRUST REGION NOT PROVIDED BY USER,
C       COMPUTE INITIAL TRUST REGION.
C
        IF(DLT.NE. (-1.D0)) GO TO 100
C       IF(DLT.EQ. (-1.))
C       THEN
          ALPHA=0.D0
          DO 10 I=1,N
            ALPHA=ALPHA+(G(I)*G(I))/(SX(I)*SX(I))
   10     CONTINUE
          BETA=0.0D0
          DO 30 I=1,N
            TMP=0.D0
            DO 20 J=I,N
              TMP=TMP + (A(J,I)*G(J))/(SX(J)*SX(J))
   20       CONTINUE
            BETA=BETA+TMP*TMP
   30     CONTINUE
          DLT=ALPHA*SQRT(ALPHA)/BETA
          DLT = MIN(DLT, STEPMX)
C$        WRITE(IPR,950)
C$        WRITE(IPR,951) ALPHA,BETA,DLT
C       ENDIF
C     ENDIF
C
  100 CONTINUE
C
C FIND NEW STEP BY MORE-HEBDON ALGORITHM
      CALL HOKSTD(NR,N,G,A,UDIAG,P,SX,RNWTLN,DLT,AMU,
     +     DLTP,PHI,PHIP0,FSTIME,SC,NWTAKE,WRK0,EPSM,IPR)
      DLTP=DLT
C
C CHECK NEW POINT AND UPDATE TRUST REGION
      CALL TRGUPD(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,
     +         DLT,IRETCD,XPLSP,FPLSP,XPLS,FPLS,MXTAKE,IPR,3,UDIAG)
      IF(IRETCD.LE.1) RETURN
      GO TO 100
C
  950 FORMAT(43H HOKDRD    INITIAL TRUST REGION NOT GIVEN. ,
     +       21H COMPUTE CAUCHY STEP.)
  951 FORMAT(18H HOKDRD    ALPHA =,E20.13/
     +       18H HOKDRD    BETA  =,E20.13/
     +       18H HOKDRD    DLT   =,E20.13)
  952 FORMAT(28H HOKDRD    CURRENT STEP (SC))
  954 FORMAT(18H0HOKDRD    RNWTLN=,E20.13)
  955 FORMAT(14H HOKDRD       ,5(E20.13,3X))
      END
      SUBROUTINE HOKSTD(NR,N,G,A,UDIAG,P,SX,RNWTLN,DLT,AMU,
     +     DLTP,PHI,PHIP0,FSTIME,SC,NWTAKE,WRK0,EPSM,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C FIND NEW STEP BY MORE-HEBDON ALGORITHM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C G(N)         --> GRADIENT AT CURRENT ITERATE, G(X)
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C P(N)         --> NEWTON STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR N
C RNWTLN       --> NEWTON STEP LENGTH
C DLT         <--> TRUST REGION RADIUS
C AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C DLTP         --> TRUST REGION RADIUS AT LAST EXIT FROM THIS ROUTINE
C PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C FSTIME      <--> BOOLEAN. =.TRUE. IF FIRST ENTRY TO THIS ROUTINE
C                  DURING K-TH ITERATION
C SC(N)       <--  CURRENT STEP
C NWTAKE      <--  BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C WRK0         --> WORKSPACE
C EPSM         --> MACHINE EPSILON
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION G(N),P(N),SX(N),SC(N),WRK0(N)
      DIMENSION A(NR,1),UDIAG(N)
      LOGICAL NWTAKE,DONE
      LOGICAL FSTIME
C
C HI AND ALO ARE CONSTANTS USED IN THIS ROUTINE.
C CHANGE HERE IF OTHER VALUES ARE TO BE SUBSTITUTED.
      IPR=IPR
      HI=1.5D0
      ALO=.75D0
C -----
      IF(RNWTLN.GT.HI*DLT) GO TO 15
C     IF(RNWTLN.LE.HI*DLT)
C     THEN
C
C       TAKE NEWTON STEP
C
        NWTAKE=.TRUE.
        DO 10 I=1,N
          SC(I)=P(I)
   10   CONTINUE
        DLT=MIN(DLT,RNWTLN)
        AMU=0.D0
C$      WRITE(IPR,951)
        RETURN
C     ELSE
C
C       NEWTON STEP NOT TAKEN
C
   15   CONTINUE
C$      WRITE(IPR,952)
        NWTAKE=.FALSE.
        IF(AMU.LE.0.D0) GO TO 20
C       IF(AMU.GT.0.)
C       THEN
          AMU=AMU- (PHI+DLTP) *((DLTP-DLT)+PHI)/(DLT*PHIP)
C$        WRITE(IPR,956) AMU
C       ENDIF
   20   CONTINUE
        PHI=RNWTLN-DLT
        IF(.NOT.FSTIME) GO TO 28
C       IF(FSTIME)
C       THEN
          DO 25 I=1,N
            WRK0(I)=SX(I)*SX(I)*P(I)
   25     CONTINUE
C
C         SOLVE L*Y = (SX**2)*P
C
          CALL FORSLD(NR,N,A,WRK0,WRK0)
          PHIP0=-DNRM2(N,WRK0,1)**2/RNWTLN
          FSTIME=.FALSE.
C       ENDIF
   28   PHIP=PHIP0
        AMULO=-PHI/PHIP
        AMUUP=0.0D0
        DO 30 I=1,N
          AMUUP=AMUUP+(G(I)*G(I))/(SX(I)*SX(I))
   30   CONTINUE
        AMUUP=SQRT(AMUUP)/DLT
        DONE=.FALSE.
C$      WRITE(IPR,956) AMU
C$      WRITE(IPR,959) PHI
C$      WRITE(IPR,960) PHIP
C$      WRITE(IPR,957) AMULO
C$      WRITE(IPR,958) AMUUP
C
C       TEST VALUE OF AMU; GENERATE NEXT AMU IF NECESSARY
C
  100   CONTINUE
        IF(DONE) RETURN
C$      WRITE(IPR,962)
        IF(AMU.GE.AMULO .AND. AMU.LE.AMUUP) GO TO 110
C       IF(AMU.LT.AMULO .OR.  AMU.GT.AMUUP)
C       THEN
          AMU=MAX(SQRT(AMULO*AMUUP),AMUUP*1.0D-3)
C$        WRITE(IPR,956) AMU
C       ENDIF
  110   CONTINUE
C
C       COPY (H,UDIAG) TO L
C       WHERE H <-- H+AMU*(SX**2) [DO NOT ACTUALLY CHANGE (H,UDIAG)]
        DO 130 J=1,N
          A(J,J)=UDIAG(J) + AMU*SX(J)*SX(J)
          IF(J.EQ.N) GO TO 130
          JP1=J+1
          DO 120 I=JP1,N
            A(I,J)=A(J,I)
  120     CONTINUE
  130   CONTINUE
C
C       FACTOR H=L(L+)
C
        CALL CHLDCD(NR,N,A,0.0D0,SQRT(EPSM),ADDMAX)
C
C       SOLVE H*P = L(L+)*SC = -G
C
        DO 140 I=1,N
          WRK0(I)=-G(I)
  140   CONTINUE
        CALL LLTSLD(NR,N,A,SC,WRK0)
C$      WRITE(IPR,955)
C$      WRITE(IPR,963) (SC(I),I=1,N)
C
C       RESET H.  NOTE SINCE UDIAG HAS NOT BEEN DESTROYED WE NEED DO
C       NOTHING HERE.  H IS IN THE UPPER PART AND IN UDIAG, STILL INTACT
C
        STEPLN=0.D0
        DO 150 I=1,N
          STEPLN=STEPLN + SX(I)*SX(I)*SC(I)*SC(I)
  150   CONTINUE
        STEPLN=SQRT(STEPLN)
        PHI=STEPLN-DLT
        DO 160 I=1,N
          WRK0(I)=SX(I)*SX(I)*SC(I)
  160   CONTINUE
        CALL FORSLD(NR,N,A,WRK0,WRK0)
        PHIP=-DNRM2(N,WRK0,1)**2/STEPLN
C$      WRITE(IPR,961) DLT,STEPLN
C$      WRITE(IPR,959) PHI
C$      WRITE(IPR,960) PHIP
        IF((ALO*DLT.GT.STEPLN .OR. STEPLN.GT.HI*DLT) .AND.
     +       (AMUUP-AMULO.GT.0.D0)) GO TO 170
C       IF((ALO*DLT.LE.STEPLN .AND. STEPLN.LE.HI*DLT) .OR.
C            (AMUUP-AMULO.LE.0.))
C       THEN
C
C         SC IS ACCEPTABLE HOKSTDEP
C
C$        WRITE(IPR,954)
          DONE=.TRUE.
          GO TO 100
C       ELSE
C
C         SC NOT ACCEPTABLE HOKSTDEP.  SELECT NEW AMU
C
  170     CONTINUE
C$        WRITE(IPR,953)
          AMULO=MAX(AMULO,AMU-(PHI/PHIP))
          IF(PHI.LT.0.D0) AMUUP=MIN(AMUUP,AMU)
          AMU=AMU-(STEPLN*PHI)/(DLT*PHIP)
C$        WRITE(IPR,956) AMU
C$        WRITE(IPR,957) AMULO
C$        WRITE(IPR,958) AMUUP
          GO TO 100
C       ENDIF
C     ENDIF
C
  951 FORMAT(27H0HOKSTD    TAKE NEWTON STEP)
  952 FORMAT(32H0HOKSTD    NEWTON STEP NOT TAKEN)
  953 FORMAT(31H HOKSTD    SC IS NOT ACCEPTABLE)
  954 FORMAT(27H HOKSTD    SC IS ACCEPTABLE)
  955 FORMAT(28H HOKSTD    CURRENT STEP (SC))
  956 FORMAT(18H HOKSTD    AMU   =,E20.13)
  957 FORMAT(18H HOKSTD    AMULO =,E20.13)
  958 FORMAT(18H HOKSTD    AMUUP =,E20.13)
  959 FORMAT(18H HOKSTD    PHI   =,E20.13)
  960 FORMAT(18H HOKSTD    PHIP  =,E20.13)
  961 FORMAT(18H HOKSTD    DLT   =,E20.13/
     +       18H HOKSTD    STEPLN=,E20.13)
  962 FORMAT(23H0HOKSTD    FIND NEW AMU)
  963 FORMAT(14H HOKSTD       ,5(E20.13,3X))
      END
      SUBROUTINE HSCHKD(NR,N,X,FCN,D1FCND,D2FCND,F,G,A,TYPSIZ,SX,RNF,
     +     ANALTL,IAGFLG,UDIAG,WRK1,WRK2,MSG,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C CHECK ANALYTIC HESSIAN AGAINST ESTIMATED HESSIAN
C  (THIS MAY BE DONE ONLY IF THE USER SUPPLIED ANALYTIC HESSIAN
C   D2FCND FILLS ONLY THE LOWER TRIANGULAR PART AND DIAGONAL OF A)
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ESTIMATE TO A ROOT OF FCN
C FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C                       FCN:  R(N) --> R(1)
C D1FCND       --> NAME OF SUBROUTINE TO EVALUATE GRADIENT OF FCN.
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C D2FCND       --> NAME OF SUBROUTINE TO EVALUATE HESSIAN OF FCN.
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C F            --> FUNCTION VALUE:  FCN(X)
C G(N)        <--  GRADIENT:  G(X)
C A(N,N)      <--  ON EXIT:  HESSIAN IN LOWER TRIANGULAR PART AND DIAG
C TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X
C SX(N)        --> DIAGONAL SCALING MATRIX:  SX(I)=1./TYPSIZ(I)
C RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN
C ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND
C                  ANALYTICAL GRADIENTS
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED
C UDIAG(N)     --> WORKSPACE
C WRK1(N)      --> WORKSPACE
C WRK2(N)      --> WORKSPACE
C MSG         <--> MESSAGE OR ERROR CODE
C                    ON INPUT : IF =1XX DO NOT COMPARE ANAL + EST HESS
C                    ON OUTPUT: =-22, PROBABLE CODING ERROR OF HESSIAN
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),G(N),A(NR,1)
      DIMENSION TYPSIZ(N),SX(N)
      DIMENSION UDIAG(N),WRK1(N),WRK2(N)
      EXTERNAL FCN,D1FCND
C
C COMPUTE FINITE DIFFERENCE APPROXIMATION A TO THE HESSIAN.
C
      IF(IAGFLG.EQ.1) CALL FSTFDD(NR,N,N,X,D1FCND,G,A,SX,RNF,WRK1,3)
      IF(IAGFLG.NE.1) CALL SNDFDD(NR,N,X,FCN,F,A,SX,RNF,WRK1,WRK2)
C
      KER=0
C
C COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART
C AND DIAGONAL OF "A" TO UDIAG
C
      DO 30 J=1,N
        UDIAG(J)=A(J,J)
        IF(J.EQ.N) GO TO 30
        JP1=J+1
        DO 25 I=JP1,N
          A(J,I)=A(I,J)
   25   CONTINUE
   30 CONTINUE
C
C COMPUTE ANALYTIC HESSIAN AND COMPARE TO FINITE DIFFERENCE
C APPROXIMATION.
C
      CALL D2FCND(NR,N,X,A)
      DO 40 J=1,N
        HS=MAX(ABS(G(J)),1.0D0)/MAX(ABS(X(J)),TYPSIZ(J))
        IF(ABS(A(J,J)-UDIAG(J)).GT.MAX(ABS(UDIAG(J)),HS)*ANALTL)
     +       KER=1
        IF(J.EQ.N) GO TO 40
        JP1=J+1
        DO 35 I=JP1,N
          IF(ABS(A(I,J)-A(J,I)).GT.MAX(ABS(A(I,J)),HS)*ANALTL) KER=1
   35   CONTINUE
   40 CONTINUE
C
      IF(KER.EQ.0) GO TO 90
        WRITE(IPR,901)
        DO 50 I=1,N
          IF(I.EQ.1) GO TO 45
          IM1=I-1
          DO 43 J=1,IM1
            WRITE(IPR,902) I,J,A(I,J),A(J,I)
   43     CONTINUE
   45     WRITE(IPR,902) I,I,A(I,I),UDIAG(I)
   50   CONTINUE
        MSG=-22
C     ENDIF
   90 CONTINUE
      RETURN
  901 FORMAT(47H HSCHKD    PROBABLE ERROR IN CODING OF ANALYTIC,
     +       18H HESSIAN FUNCTION./
     +       21H HSCHKD      ROW  COL,14X,8HANALYTIC,14X,10H(ESTIMATE))
  902 FORMAT(11H HSCHKD    ,2I5,2X,E20.13,2X,1H(,E20.13,1H))
      END
      SUBROUTINE HSNNTD(NR,N,A,SX,METHOD)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C PROVIDE INITIAL HESSIAN WHEN USING SECANT UPDATES
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--  INITIAL HESSIAN (LOWER TRIANGULAR MATRIX)
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1,2 FACTORED SECANT METHOD USED
C                    =3   UNFACTORED SECANT METHOD USED
C
      DIMENSION A(NR,1),SX(N)
C
      DO 100 J=1,N
        IF(METHOD.EQ.3) A(J,J)=SX(J)*SX(J)
        IF(METHOD.NE.3) A(J,J)=SX(J)
        IF(J.EQ.N) GO TO 100
        JP1=J+1
        DO 90 I=JP1,N
          A(I,J)=0.D0
   90   CONTINUE
  100 CONTINUE
      RETURN
      END
      INTEGER FUNCTION I1MACH(I)

c*********************************************************************72
c
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  840405   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS) 
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C   These machine constant routines must be activated for
C   a particular environment.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function 
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for 
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers. 
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude. 
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL 
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
C
C      DATA IMACH( 1) /    5 / 
C      DATA IMACH( 2) /    6 / 
C      DATA IMACH( 3) /    7 / 
C      DATA IMACH( 4) /    6 / 
C      DATA IMACH( 5) /   60 / 
C      DATA IMACH( 6) /   10 / 
C      DATA IMACH( 7) /    2 / 
C      DATA IMACH( 8) /   48 / 
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 / 
C      DATA IMACH(11) /   48 / 
C      DATA IMACH(12) / -974 / 
C      DATA IMACH(13) / 1070 / 
C      DATA IMACH(14) /   96 / 
C      DATA IMACH(15) / -927 / 
C      DATA IMACH(16) / 1070 / 
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
C
C     DATA IMACH( 1) /      5 /
C     DATA IMACH( 2) /      6 /
C     DATA IMACH( 3) /      7 /
C     DATA IMACH( 4) /      6 /
C     DATA IMACH( 5) /     64 /
C     DATA IMACH( 6) /      8 /
C     DATA IMACH( 7) /      2 /
C     DATA IMACH( 8) /     47 /
C     DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
C     DATA IMACH(10) /      2 /
C     DATA IMACH(11) /     47 /
C     DATA IMACH(12) / -28625 /
C     DATA IMACH(13) /  28718 /
C     DATA IMACH(14) /     94 /
C     DATA IMACH(15) / -28625 /
C     DATA IMACH(16) /  28718 /
C
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA IMACH( 1) /    5 / 
C     DATA IMACH( 2) /    6 / 
C     DATA IMACH( 3) /    7 / 
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 / 
C     DATA IMACH( 6) /   10 / 
C     DATA IMACH( 7) /    2 / 
C     DATA IMACH( 8) /   48 / 
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 / 
C     DATA IMACH(11) /   47 / 
C     DATA IMACH(12) / -929 / 
C     DATA IMACH(13) / 1070 / 
C     DATA IMACH(14) /   94 / 
C     DATA IMACH(15) / -929 / 
C     DATA IMACH(16) / 1069 / 
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
C
      DATA IMACH/5,6,0,6,32,4,2,31,2147483647,2,24,
     * -125,127,53,-1021,1023/
C               NOTE! I1MACH(3) IS NOT WELL DEFINED AND IS SET TO ZERO.
C
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    5 / 
C     DATA IMACH( 2) /    6 / 
C     DATA IMACH( 3) /    5 / 
C     DATA IMACH( 4) /    6 / 
C     DATA IMACH( 5) /   36 / 
C     DATA IMACH( 6) /    5 / 
C     DATA IMACH( 7) /    2 / 
C     DATA IMACH( 8) /   35 / 
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 / 
C     DATA IMACH(11) /   27 / 
C     DATA IMACH(12) / -128 / 
C     DATA IMACH(13) /  127 / 
C     DATA IMACH(14) /   54 / 
C     DATA IMACH(15) / -101 / 
C     DATA IMACH(16) /  127 / 
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    5 / 
C     DATA IMACH( 2) /    6 / 
C     DATA IMACH( 3) /    5 / 
C     DATA IMACH( 4) /    6 / 
C     DATA IMACH( 5) /   36 / 
C     DATA IMACH( 6) /    5 / 
C     DATA IMACH( 7) /    2 / 
C     DATA IMACH( 8) /   35 / 
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 / 
C     DATA IMACH(11) /   27 / 
C     DATA IMACH(12) / -128 / 
C     DATA IMACH(13) /  127 / 
C     DATA IMACH(14) /   62 / 
C     DATA IMACH(15) / -128 / 
C     DATA IMACH(16) /  127 / 
C
C
C     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
C       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
C       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    6 /
C      DATA IMACH( 4) /    0 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16) 
     1   CALL XERROR ( 'I1MACH -- I OUT OF BOUNDS',25,1,2)
C
      I1MACH=IMACH(I)
      RETURN
C
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  D1A2
C***KEYWORDS  LIBRARY=SLATEC(BLAS),
C             TYPE=DOUBLE PRECISION(ISAMAX-S IDAMAX-D ICAMAX-C),
C             LINEAR ALGEBRA,MAXIMUM COMPONENT,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find the smallest index of that component of a d.p. vector
C            having the maximum magnitude.
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  IDAMAX
C
      DOUBLE PRECISION DX(*),DMAX,XMAG
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAX = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END
      FUNCTION INITDS(DOS,NOS,ETA)

c*********************************************************************72
c
C***BEGIN PROLOGUE  INITDS
C***DATE WRITTEN   770601   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  C3A2
C***KEYWORDS  CHEBYSHEV,DOUBLE PRECISION,INITIALIZE,
C             ORTHOGONAL POLYNOMIAL,SERIES,SPECIAL FUNCTION
C***AUTHOR  FULLERTON, W., (LANL)
C***PURPOSE  Initializes the d.p. properly normalized orthogonal
C            polynomial series to determine the number of terms needed
C            for specific accuracy.
C***DESCRIPTION
C
C Initialize the double precision orthogonal series DOS so that INITDS
C is the number of terms needed to insure the error is no larger than
C ETA.  Ordinarily ETA will be chosen to be one-tenth machine precision
C
C             Input Arguments --
C DOS    dble prec array of NOS coefficients in an orthogonal series.
C NOS    number of coefficients in DOS.
C ETA    requested accuracy of series.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  INITDS
C
      DOUBLE PRECISION DOS(NOS)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS.LT.1) CALL XERROR ( 'INITDS  NUMBER OF COEFFICIENTS LT 1',
     1 35, 2, 2)
C
      ERR = 0.
      DO 10 II=1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(SNGL(DOS(I)))
        IF (ERR.GT.ETA) GO TO 20
 10   CONTINUE
C
 20   IF (I.EQ.NOS) CALL XERROR ( 'INITDS  ETA MAY BE TOO SMALL', 28,
     1  1, 2)
      INITDS = I
C
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)

c*********************************************************************72
c
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE LLTSLD(NR,N,A,X,B)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C SOLVE AX=B WHERE A HAS THE FORM L(L-TRANSPOSE)
C BUT ONLY THE LOWER TRIANGULAR PART, L, IS STORED.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> MATRIX OF FORM L(L-TRANSPOSE).
C                  ON RETURN A IS UNCHANGED.
C X(N)        <--  SOLUTION VECTOR
C B(N)         --> RIGHT-HAND SIDE VECTOR
C
C NOTE
C ----
C IF B IS NOT REQUIRED BY CALLING PROGRAM, THEN
C B AND X MAY SHARE THE SAME STORAGE.
C
      DIMENSION A(NR,1),X(N),B(N)
C
C FORWARD SOLVE, RESULT IN X
C
      CALL FORSLD(NR,N,A,X,B)
C
C BACK SOLVE, RESULT IN X
C
      CALL BAKSLD(NR,N,A,X,X)
      RETURN
      END
      SUBROUTINE LNSRCD(N,X,F,G,P,XPLS,FPLS,FCN,MXTAKE,
     +   IRETCD,STEPMX,STEPTL,SX,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PURPOSE
C -------
C FIND A NEXT NEWTON ITERATE BY LINE SEARCH.
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE:   X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C P(N)         --> NON-ZERO NEWTON STEP
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C IRETCD      <--  RETURN CODE
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
C INTERNAL VARIABLES
C ------------------
C SLN              NEWTON LENGTH
C RLN              RELATIVE LENGTH OF NEWTON STEP
C
      INTEGER N,IRETCD
      DIMENSION SX(N)
      DIMENSION X(N),G(N),P(N)
      DIMENSION XPLS(N)
      LOGICAL MXTAKE
C
      IPR=IPR
      MXTAKE=.FALSE.
      IRETCD=2
C$    WRITE(IPR,954)
C$    WRITE(IPR,955) (P(I),I=1,N)
      TMP=0.0D0
      DO 5 I=1,N
        TMP=TMP+SX(I)*SX(I)*P(I)*P(I)
    5 CONTINUE
      SLN=SQRT(TMP)
      IF(SLN.LE.STEPMX) GO TO 10
C
C NEWTON STEP LONGER THAN MAXIMUM ALLOWED
        SCL=STEPMX/SLN
        CALL SCLMLD(N,SCL,P,P)
        SLN=STEPMX
C$      WRITE(IPR,954)
C$      WRITE(IPR,955) (P(I),I=1,N)
   10 CONTINUE
      SLP=DDOT(N,G,1,P,1)
      RLN=0.D0
      DO 15 I=1,N
        RLN=MAX(RLN,ABS(P(I))/MAX(ABS(X(I)),1.D0/SX(I)))
   15 CONTINUE
      RMNLMB=STEPTL/RLN
      ALMBDA=1.0D0
C$    WRITE(IPR,952) SLN,SLP,RMNLMB,STEPMX,STEPTL
C
C LOOP
C CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY.
C
  100 CONTINUE
      IF(IRETCD.LT.2) RETURN
      DO 105 I=1,N
        XPLS(I)=X(I) + ALMBDA*P(I)
  105 CONTINUE
      CALL FCN(N,XPLS,FPLS)
C$    WRITE(IPR,950) ALMBDA
C$    WRITE(IPR,951)
C$    WRITE(IPR,955) (XPLS(I),I=1,N)
C$    WRITE(IPR,953) FPLS
      IF(FPLS.GT. F+SLP*1.D-4*ALMBDA) GO TO 130
C     IF(FPLS.LE. F+SLP*1.D-4*ALMBDA)
C     THEN
C
C SOLUTION FOUND
C
        IRETCD=0
        IF(ALMBDA.EQ.1.0D0 .AND. SLN.GT. .99D0*STEPMX) MXTAKE=.TRUE.
        GO TO 100
C
C SOLUTION NOT (YET) FOUND
C
C     ELSE
  130   IF(ALMBDA .GE. RMNLMB) GO TO 140
C       IF(ALMBDA .LT. RMNLMB)
C       THEN
C
C NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X
C
          IRETCD=1
          GO TO 100
C       ELSE
C
C CALCULATE NEW LAMBDA
C
  140     IF(ALMBDA.NE.1.0D0) GO TO 150
C         IF(ALMBDA.EQ.1.0)
C         THEN
C
C FIRST BACKTRACK: QUADRATIC FIT
C
            TLMBDA=-SLP/(2.D0*(FPLS-F-SLP))
            GO TO 170
C         ELSE
C
C ALL SUBSEQUENT BACKTRACKS: CUBIC FIT
C
  150       T1=FPLS-F-ALMBDA*SLP
            T2=PFPLS-F-PLMBDA*SLP
            T3=1.0D0/(ALMBDA-PLMBDA)
            A=T3*(T1/(ALMBDA*ALMBDA) - T2/(PLMBDA*PLMBDA))
            B=T3*(T2*ALMBDA/(PLMBDA*PLMBDA)
     +           - T1*PLMBDA/(ALMBDA*ALMBDA) )
            DISC=B*B-3.0D0*A*SLP
            IF(DISC.LE. B*B) GO TO 160
C           IF(DISC.GT. B*B)
C           THEN
C
C ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM
C
              TLMBDA=(-B+SIGN(1.0D0,A)*SQRT(DISC))/(3.0D0*A)
              GO TO 165
C           ELSE
C
C BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM
C
  160         TLMBDA=(-B-SIGN(1.0D0,A)*SQRT(DISC))/(3.0D0*A)
C           ENDIF
  165       IF(TLMBDA.GT. .5D0*ALMBDA) TLMBDA=.5D0*ALMBDA
C         ENDIF
  170     PLMBDA=ALMBDA
          PFPLS=FPLS
          IF(TLMBDA.GE. ALMBDA*.1D0) GO TO 180
C         IF(TLMBDA.LT.ALMBDA/10.)
C         THEN
            ALMBDA=ALMBDA*.1D0
            GO TO 190
C         ELSE
  180       ALMBDA=TLMBDA
C         ENDIF
C       ENDIF
C     ENDIF
  190 GO TO 100
  950 FORMAT(18H LNSRCD    ALMBDA=,E20.13)
  951 FORMAT(29H LNSRCD    NEW ITERATE (XPLS))
  952 FORMAT(18H LNSRCD    SLN   =,E20.13/
     +       18H LNSRCD    SLP   =,E20.13/
     +       18H LNSRCD    RMNLMB=,E20.13/
     +       18H LNSRCD    STEPMX=,E20.13/
     +       18H LNSRCD    STEPTL=,E20.13)
  953 FORMAT(19H LNSRCD    F(XPLS)=,E20.13)
  954 FORMAT(26H0LNSRCD    NEWTON STEP (P))
  955 FORMAT(14H LNSRCD       ,5(E20.13,3X))
      END
      SUBROUTINE MVMLLD(NR,N,A,X,Y)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C COMPUTE Y=LX
C WHERE L IS A LOWER TRIANGULAR MATRIX STORED IN A
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> LOWER TRIANGULAR (N*N) MATRIX
C X(N)         --> OPERAND VECTOR
C Y(N)        <--  RESULT VECTOR
C
C NOTE
C ----
C X AND Y CANNOT SHARE STORAGE
C
      DIMENSION A(NR,1),X(N),Y(N)
      DO 30 I=1,N
        SUM=0.D0
        DO 10 J=1,I
          SUM=SUM+A(I,J)*X(J)
   10   CONTINUE
        Y(I)=SUM
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MVMLSD(NR,N,A,X,Y)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C COMPUTE Y=AX
C WHERE "A" IS A SYMMETRIC (N*N) MATRIX STORED IN ITS LOWER
C TRIANGULAR PART AND X,Y ARE N-VECTORS
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)       --> SYMMETRIC (N*N) MATRIX STORED IN
C                  LOWER TRIANGULAR PART AND DIAGONAL
C X(N)         --> OPERAND VECTOR
C Y(N)        <--  RESULT VECTOR
C
C NOTE
C ----
C X AND Y CANNOT SHARE STORAGE.
C
      DIMENSION A(NR,1),X(N),Y(N)
      DO 30 I=1,N
        SUM=0.D0
        DO 10 J=1,I
          SUM=SUM+A(I,J)*X(J)
   10   CONTINUE
        IF(I.EQ.N) GO TO 25
        IP1=I+1
        DO 20 J=IP1,N
          SUM=SUM+A(J,I)*X(J)
   20   CONTINUE
   25   Y(I)=SUM
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MVMLUD(NR,N,A,X,Y)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C COMPUTE Y=(L+)X
C WHERE L IS A LOWER TRIANGULAR MATRIX STORED IN A
C (L-TRANSPOSE (L+) IS TAKEN IMPLICITLY)
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(NR,1)       --> LOWER TRIANGULAR (N*N) MATRIX
C X(N)         --> OPERAND VECTOR
C Y(N)        <--  RESULT VECTOR
C
C NOTE
C ----
C X AND Y CANNOT SHARE STORAGE
C
      DIMENSION A(NR,1),X(N),Y(N)
      DO 30 I=1,N
        SUM=0.D0
        DO 10 J=I,N
          SUM=SUM+A(J,I)*X(J)
   10   CONTINUE
        Y(I)=SUM
   30 CONTINUE
      RETURN
      END
      FUNCTION NUMXER(NERR)

c*********************************************************************72
c
C***BEGIN PROLOGUE  NUMXER
C***REFER TO  XERROR
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        NUMXER returns the most recent error number,
C        in both NUMXER and the parameter NERR.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 JUNE 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  NUMXER
C***FIRST EXECUTABLE STATEMENT  NUMXER
      NERR = J4SAVE(1,0,.FALSE.)
      NUMXER = NERR
      RETURN
      END
      SUBROUTINE OPTCHD(N,X,TYPSIZ,SX,FSCALE,GRADTL,ITNLIM,NDIGIT,EPSM,
     +     DLT,METHOD,IEXP,IAGFLG,IAHFLG,STEPMX,MSG,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C CHECK INPUT FOR REASONABLENESS
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ON ENTRY, ESTIMATE TO ROOT OF FCN
C TYPSIZ(N)   <--> TYPICAL SIZE OF EACH COMPONENT OF X
C SX(N)       <--  DIAGONAL SCALING MATRIX FOR X
C FSCALE      <--> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN
C GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C ITNLIM      <--> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C NDIGIT      <--> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN
C EPSM         --> MACHINE EPSILON
C DLT         <--> TRUST REGION RADIUS
C METHOD      <--> ALGORITHM INDICATOR
C IEXP        <--> EXPENSE FLAG
C IAGFLG      <--> =1 IF ANALYTIC GRADIENT SUPPLIED
C IAHFLG      <--> =1 IF ANALYTIC HESSIAN SUPPLIED
C STEPMX      <--> MAXIMUM STEP SIZE
C MSG         <--> MESSAGE AND ERROR CODE
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),TYPSIZ(N),SX(N)
C
C CHECK THAT PARAMETERS ONLY TAKE ON ACCEPTABLE VALUES.
C IF NOT, SET THEM TO DEFAULT VALUES.
      IF(METHOD.LT.1 .OR. METHOD.GT.3) METHOD=1
      IF(IAGFLG.NE.1) IAGFLG=0
      IF(IAHFLG.NE.1) IAHFLG=0
      IF(IEXP.NE.0) IEXP=1
      IF(MOD(MSG/2,2).EQ.1 .AND. IAGFLG.EQ.0) GO TO 830
      IF(MOD(MSG/4,2).EQ.1 .AND. IAHFLG.EQ.0) GO TO 835
C
C CHECK DIMENSION OF PROBLEM
C
      IF(N.LE.0) GO TO 805
      IF(N.EQ.1 .AND. MOD(MSG,2).EQ.0) GO TO 810
C
C COMPUTE SCALE MATRIX
C
      DO 10 I=1,N
        IF(TYPSIZ(I).EQ.0.D0) TYPSIZ(I)=1.0D0
        IF(TYPSIZ(I).LT.0.D0) TYPSIZ(I)=-TYPSIZ(I)
        SX(I)=1.0D0/TYPSIZ(I)
   10 CONTINUE
C
C CHECK MAXIMUM STEP SIZE
C
      IF (STEPMX .GT. 0.0D0) GO TO 20
      STPSIZ = 0.0D0
      DO 15 I = 1, N
         STPSIZ = STPSIZ + X(I)*X(I)*SX(I)*SX(I)
   15 CONTINUE
      STPSIZ = SQRT(STPSIZ)
      STEPMX = MAX(1.0D3*STPSIZ, 1.0D3)
   20 CONTINUE
C CHECK FUNCTION SCALE
      IF(FSCALE.EQ.0.D0) FSCALE=1.0D0
      IF(FSCALE.LT.0.D0) FSCALE=-FSCALE
C
C CHECK GRADIENT TOLERANCE
      IF(GRADTL.LT.0.D0) GO TO 815
C
C CHECK ITERATION LIMIT
      IF(ITNLIM.LE.0) GO TO 820
C
C CHECK NUMBER OF DIGITS OF ACCURACY IN FUNCTION FCN
      IF(NDIGIT.EQ.0) GO TO 825
      IF(NDIGIT.LT.0) NDIGIT=-LOG10(EPSM)
C
C CHECK TRUST REGION RADIUS
      IF(DLT.LE.0.D0) DLT=-1.0D0
      IF (DLT .GT. STEPMX) DLT = STEPMX
      RETURN
C
C ERROR EXITS
C
  805 WRITE(IPR,901) N
      MSG=-1
      GO TO 895
  810 WRITE(IPR,902)
      MSG=-2
      GO TO 895
  815 WRITE(IPR,903) GRADTL
      MSG=-3
      GO TO 895
  820 WRITE(IPR,904) ITNLIM
      MSG=-4
      GO TO 895
  825 WRITE(IPR,905) NDIGIT
      MSG=-5
      GO TO 895
  830 WRITE(IPR,906) MSG,IAGFLG
      MSG=-6
      GO TO 895
  835 WRITE(IPR,907) MSG,IAHFLG
      MSG=-7
  895 RETURN
  901 FORMAT(32H0OPTCHD    ILLEGAL DIMENSION, N=,I5)
  902 FORMAT(55H0OPTCHD    +++ WARNING +++  THIS PACKAGE IS INEFFICIENT,
     +       26H FOR PROBLEMS OF SIZE N=1./
     +       48H OPTCHD    CHECK INSTALLATION LIBRARIES FOR MORE,
     +       22H APPROPRIATE ROUTINES./
     +       41H OPTCHD    IF NONE, SET MSG AND RESUBMIT.)
  903 FORMAT(38H0OPTCHD    ILLEGAL TOLERANCE.  GRADTL=,E20.13)
  904 FORMAT(44H0OPTCHD    ILLEGAL ITERATION LIMIT.  ITNLIM=,I5)
  905 FORMAT(52H0OPTCHD    MINIMIZATION FUNCTION HAS NO GOOD DIGITS.,
     +        9H  NDIGIT=,I5)
  906 FORMAT(50H0OPTCHD    USER REQUESTS THAT ANALYTIC GRADIENT BE,
     +       33H ACCEPTED AS PROPERLY CODED (MSG=,I5, 2H),/
     +       45H OPTCHD    BUT ANALYTIC GRADIENT NOT SUPPLIED,
     +        9H (IAGFLG=,I5, 2H).)
  907 FORMAT(49H0OPTCHD    USER REQUESTS THAT ANALYTIC HESSIAN BE,
     +       33H ACCEPTED AS PROPERLY CODED (MSG=,I5, 2H),/
     +       44H OPTCHD    BUT ANALYTIC HESSIAN NOT SUPPLIED,
     +        9H (IAHFLG=,I5, 2H).)
      END
      SUBROUTINE OPTDRD(NR,N,X,FCN,D1FCND,D2FCND,TYPSIZ,FSCALE,
     +     METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
     +     DLT,GRADTL,STEPMX,STEPTL,
     +     XPLS,FPLS,GPLS,ITRMCD,
     +     A,UDIAG,G,P,SX,WRK0,WRK1,WRK2,WRK3)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C DRIVER FOR NON-LINEAR OPTIMIZATION PROBLEM
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ON ENTRY: ESTIMATE TO A ROOT OF FCN
C FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C                            FCN: R(N) --> R(1)
C D1FCND       --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE GRADIENT
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C D2FCND       --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE HESSIAN OF
C                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE
C TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1 LINE SEARCH
C                    =2 DOUBLE DOGLEG
C                    =3 MORE-HEBDON
C IEXP         --> =1 IF OPTIMIZATION FUNCTION FCN IS EXPENSIVE TO
C                  EVALUATE, =0 OTHERWISE.  IF SET THEN HESSIAN WILL
C                  BE EVALUATED BY SECANT UPDATE INSTEAD OF
C                  ANALYTICALLY OR BY FINITE DIFFERENCES
C MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN
C                    AUTOMATIC CHECKS
C                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR
C NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN
C ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED
C IAHFLG       --> =1 IF ANALYTIC HESSIAN SUPPLIED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C DLT          --> TRUST REGION RADIUS
C GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C XPLS(N)     <--> ON EXIT:  XPLS IS LOCAL MINIMUM
C FPLS        <--> ON EXIT:  FUNCTION VALUE AT SOLUTION, XPLS
C GPLS(N)     <--> ON EXIT:  GRADIENT AT SOLUTION XPLS
C ITRMCD      <--  TERMINATION CODE
C A(N,N)       --> WORKSPACE FOR HESSIAN (OR ESTIMATE)
C                  AND ITS CHOLESKY DECOMPOSITION
C UDIAG(N)     --> WORKSPACE [FOR DIAGONAL OF HESSIAN]
C G(N)         --> WORKSPACE (FOR GRADIENT AT CURRENT ITERATE)
C P(N)         --> WORKSPACE FOR STEP
C SX(N)        --> WORKSPACE (FOR DIAGONAL SCALING MATRIX)
C WRK0(N)      --> WORKSPACE
C WRK1(N)      --> WORKSPACE
C WRK2(N)      --> WORKSPACE
C WRK3(N)      --> WORKSPACE
C
C
C INTERNAL VARIABLES
C ------------------
C ANALTL           TOLERANCE FOR COMPARISON OF ESTIMATED AND
C                  ANALYTICAL GRADIENTS AND HESSIANS
C EPSM             MACHINE EPSILON
C F                FUNCTION VALUE: FCN(X)
C ITNCNT           CURRENT ITERATION, K
C RNF              RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN.
C                       NOISE=10.**(-NDIGIT)
C
      DIMENSION X(N),XPLS(N),G(N),GPLS(N),P(N)
      DIMENSION TYPSIZ(N),SX(N)
      DIMENSION A(NR,1),UDIAG(N)
      DIMENSION WRK0(N),WRK1(N),WRK2(N),WRK3(N)
      LOGICAL MXTAKE,NOUPDT
      EXTERNAL FCN,D1FCND,D2FCND
C
C INITIALIZATION
C --------------
      DO 10 I=1,N
        P(I)=0.D0
   10 CONTINUE
      ITNCNT=0
      IRETCD=-1
      EPSM=D1MACH(4)
      CALL OPTCHD(N,X,TYPSIZ,SX,FSCALE,GRADTL,ITNLIM,NDIGIT,EPSM,
     +     DLT,METHOD,IEXP,IAGFLG,IAHFLG,STEPMX,MSG,IPR)
      IF(MSG.LT.0) RETURN
      RNF=MAX(10.0D0**(-NDIGIT),EPSM)
      ANALTL=MAX(1.0D-2,SQRT(RNF))
C
      IF(MOD(MSG/8,2).EQ.1) GO TO 15
      WRITE(IPR,901)
      WRITE(IPR,900) (TYPSIZ(I),I=1,N)
      WRITE(IPR,902)
      WRITE(IPR,900) (SX(I),I=1,N)
      WRITE(IPR,903) FSCALE
      WRITE(IPR,904) NDIGIT,IAGFLG,IAHFLG,IEXP,METHOD,ITNLIM,EPSM
      WRITE(IPR,905) STEPMX,STEPTL,GRADTL,DLT,RNF,ANALTL
   15 CONTINUE
C
C EVALUATE FCN(X)
C
      CALL FCN(N,X,F)
C
C EVALUATE ANALYTIC OR FINITE DIFFERENCE GRADIENT AND CHECK ANALYTIC
C GRADIENT, IF REQUESTED.
C
      IF (IAGFLG .EQ. 1) GO TO 20
C     IF (IAGFLG .EQ. 0)
C     THEN
        CALL FSTFDD (1, 1, N, X, FCN, F, G, SX, RNF, WRK, 1)
        GO TO 25
C
   20 CALL D1FCND (N, X, G)
      IF (MOD(MSG/2,2) .EQ. 1) GO TO 25
C     IF (MOD(MSG/2,2).EQ.0)
C     THEN
        CALL GRCHKD (N, X, FCN, F, G, TYPSIZ, SX, FSCALE,
     1    RNF, ANALTL, WRK1, MSG, IPR)
        IF (MSG .LT. 0) RETURN
   25 CONTINUE
C
      CALL OPTSTD(N,X,F,G,WRK1,ITNCNT,ICSCMX,
     +            ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,MXTAKE,
     +            IPR,MSG)
      IF(ITRMCD.NE.0) GO TO 700
C
      IF(IEXP.NE.1) GO TO 80
C
C IF OPTIMIZATION FUNCTION EXPENSIVE TO EVALUATE (IEXP=1), THEN
C HESSIAN WILL BE OBTAINED BY SECANT UPDATES.  GET INITIAL HESSIAN.
C
      CALL HSNNTD(NR,N,A,SX,METHOD)
      GO TO 90
   80 CONTINUE
C
C EVALUATE ANALYTIC OR FINITE DIFFERENCE HESSIAN AND CHECK ANALYTIC
C HESSIAN IF REQUESTED (ONLY IF USER-SUPPLIED ANALYTIC HESSIAN
C ROUTINE D2FCND FILLS ONLY LOWER TRIANGULAR PART AND DIAGONAL OF A).
C
      IF (IAHFLG .EQ. 1) GO TO 82
C     IF (IAHFLG .EQ. 0)
C     THEN
         IF (IAGFLG .EQ. 1) CALL FSTFDD (NR, N, N, X, D1FCND, G, A, SX,
     1      RNF, WRK1, 3)
         IF (IAGFLG .NE. 1) CALL SNDFDD (NR, N, X, FCN, F, A, SX, RNF,
     1      WRK1, WRK2)
         GO TO 88
C
C     ELSE
   82    IF (MOD(MSG/4,2).EQ.0) GO TO 85
C        IF (MOD(MSG/4, 2) .EQ. 1)
C        THEN
            CALL D2FCND (NR, N, X, A)
            GO TO 88
C
C        ELSE
   85       CALL HSCHKD (NR, N, X, FCN, D1FCND, D2FCND, F, G, A, TYPSIZ,
     1         SX, RNF, ANALTL, IAGFLG, UDIAG, WRK1, WRK2, MSG, IPR)
C
C           HSCHKD EVALUATES D2FCND AND CHECKS IT AGAINST THE FINITE
C           DIFFERENCE HESSIAN WHICH IT CALCULATES BY CALLING FSTFDD
C           (IF IAGFLG .EQ. 1) OR SNDFDD (OTHERWISE).
C
            IF (MSG .LT. 0) RETURN
   88 CONTINUE
C
   90 IF(MOD(MSG/8,2).EQ.0)
     +     CALL RESLTD(NR,N,X,F,G,A,P,ITNCNT,1,IPR)
C
C
C ITERATION
C ---------
  100 ITNCNT=ITNCNT+1
C
C FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION
C (SKIP THIS STEP IF LINE SEARCH OR DOGSTEP TECHNIQUES BEING USED WITH
C SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALREADY OBTAINED FROM
C SECFCD.)
C
      IF(IEXP.EQ.1 .AND. METHOD.NE.3) GO TO 105
  103   CALL CHLHSD(NR,N,A,EPSM,SX,UDIAG)
  105 CONTINUE
C
C SOLVE FOR NEWTON STEP:  AP=-G
C
      DO 110 I=1,N
        WRK1(I)=-G(I)
  110 CONTINUE
      CALL LLTSLD(NR,N,A,P,WRK1)
C
C DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS=X + P
C OR TO CHOOSE XPLS BY A GLOBAL STRATEGY.
C
      IF (IAGFLG .NE. 0 .OR. METHOD .EQ. 1) GO TO 111
      DLTSAV = DLT
      IF (METHOD .EQ. 2) GO TO 111
      AMUSAV = AMU
      DLPSAV = DLTP
      PHISAV = PHI
      PHPSAV = PHIP0
  111 IF(METHOD.EQ.1)
     +     CALL LNSRCD(N,X,F,G,P,XPLS,FPLS,FCN,MXTAKE,IRETCD,
     +     STEPMX,STEPTL,SX,IPR)
      IF(METHOD.EQ.2)
     +     CALL DGDRVD(NR,N,X,F,G,A,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,WRK0,WRK1,WRK2,WRK3,IPR)
      IF(METHOD.EQ.3)
     +     CALL HOKDRD(NR,N,X,F,G,A,UDIAG,P,XPLS,FPLS,FCN,SX,STEPMX,
     +     STEPTL,DLT,IRETCD,MXTAKE,AMU,DLTP,PHI,PHIP0,WRK0,
     +     WRK1,WRK2,EPSM,ITNCNT,IPR)
C
C IF COULD NOT FIND SATISFACTORY STEP AND FORWARD DIFFERENCE
C GRADIENT WAS USED, RETRY USING CENTRAL DIFFERENCE GRADIENT.
C
      IF (IRETCD .NE. 1 .OR. IAGFLG .NE. 0) GO TO 112
C     IF (IRETCD .EQ. 1 .AND. IAGFLG .EQ. 0)
C     THEN
C
C        SET IAGFLG FOR CENTRAL DIFFERENCES
C
         IAGFLG = -1
         WRITE(*,906) ITNCNT
C
         CALL FSTCDD (N, X, FCN, SX, RNF, G)
         IF (METHOD .EQ. 1) GO TO 105
         DLT = DLTSAV
         IF (METHOD .EQ. 2) GO TO 105
         AMU = AMUSAV
         DLTP = DLPSAV
         PHI = PHISAV
         PHIP0 = PHPSAV
         GO TO 103
C     ENDIF
C
C CALCULATE STEP FOR OUTPUT
C
  112 CONTINUE
      DO 114 I = 1, N
         P(I) = XPLS(I) - X(I)
  114 CONTINUE
C
C CALCULATE GRADIENT AT XPLS
C
      IF (IAGFLG .EQ. (-1)) GO TO 116
      IF (IAGFLG .EQ. 0) GO TO 118
C
C ANALYTIC GRADIENT
      CALL D1FCND (N, XPLS, GPLS)
      GO TO 120
C
C CENTRAL DIFFERENCE GRADIENT
  116 CALL FSTCDD (N, XPLS, FCN, SX, RNF, GPLS)
      GO TO 120
C
C FORWARD DIFFERENCE GRADIENT
  118 CALL FSTFDD (1, 1, N, XPLS, FCN, FPLS, GPLS, SX, RNF, WRK, 1)
  120 CONTINUE
C
C CHECK WHETHER STOPPING CRITERIA SATISFIED
C
      CALL OPTSTD(N,XPLS,FPLS,GPLS,X,ITNCNT,ICSCMX,
     +            ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,MXTAKE,
     +            IPR,MSG)
      IF(ITRMCD.NE.0) GO TO 690
C
C EVALUATE HESSIAN AT XPLS
C
      IF(IEXP.EQ.0) GO TO 130
      IF(METHOD.EQ.3)
     +     CALL SECNFD(NR,N,X,G,A,UDIAG,XPLS,GPLS,EPSM,ITNCNT,RNF,
     +     IAGFLG,NOUPDT,WRK1,WRK2,WRK3)
      IF(METHOD.NE.3)
     +     CALL SECFCD(NR,N,X,G,A,XPLS,GPLS,EPSM,ITNCNT,RNF,IAGFLG,
     +     NOUPDT,WRK0,WRK1,WRK2,WRK3)
      GO TO 150
  130 IF(IAHFLG.EQ.1) GO TO 140
      IF(IAGFLG.EQ.1)
     +     CALL FSTFDD(NR,N,N,XPLS,D1FCND,GPLS,A,SX,RNF,WRK1,3)
      IF(IAGFLG.NE.1) CALL SNDFDD(NR,N,XPLS,FCN,FPLS,A,SX,RNF,WRK1,WRK2)
      GO TO 150
  140 CALL D2FCND(NR,N,XPLS,A)
  150 CONTINUE
      IF(MOD(MSG/16,2).EQ.1)
     +     CALL RESLTD(NR,N,XPLS,FPLS,GPLS,A,P,ITNCNT,1,IPR)
C
C X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS
C
      F=FPLS
      DO 160 I=1,N
        X(I)=XPLS(I)
        G(I)=GPLS(I)
  160 CONTINUE
      GO TO 100
C
C TERMINATION
C -----------
C RESET XPLS,FPLS,GPLS,  IF PREVIOUS ITERATE SOLUTION
C
  690 IF(ITRMCD.NE.3) GO TO 710
  700 CONTINUE
      FPLS=F
      DO 705 I=1,N
        XPLS(I)=X(I)
        GPLS(I)=G(I)
  705 CONTINUE
C
C PRINT RESULTS
C
  710 CONTINUE
      IF(MOD(MSG/8,2).EQ.0)
     +     CALL RESLTD(NR,N,XPLS,FPLS,GPLS,A,P,ITNCNT,0,IPR)
      MSG=0
      RETURN
C
  900 FORMAT(14H OPTDRD       ,5(E20.13,3X))
  901 FORMAT(20H0OPTDRD    TYPICAL X)
  902 FORMAT(40H OPTDRD    DIAGONAL SCALING MATRIX FOR X)
  903 FORMAT(22H OPTDRD    TYPICAL F =,E20.13)
  904 FORMAT(40H0OPTDRD    NUMBER OF GOOD DIGITS IN FCN=,I5/
     +       27H OPTDRD    GRADIENT FLAG  =,I5,18H   (=1 IF ANALYTIC,
     +       19H GRADIENT SUPPLIED)/
     +       27H OPTDRD    HESSIAN FLAG   =,I5,18H   (=1 IF ANALYTIC,
     +       18H HESSIAN SUPPLIED)/
     +       27H OPTDRD    EXPENSE FLAG   =,I5, 9H   (=1 IF,
     +       45H MINIMIZATION FUNCTION EXPENSIVE TO EVALUATE)/
     +       27H OPTDRD    METHOD TO USE  =,I5,19H   (=1,2,3 FOR LINE,
     +       49H SEARCH, DOUBLE DOGLEG, MORE-HEBDON RESPECTIVELY)/
     +       27H OPTDRD    ITERATION LIMIT=,I5/
     +       27H OPTDRD    MACHINE EPSILON=,E20.13)
  905 FORMAT(30H0OPTDRD    MAXIMUM STEP SIZE =,E20.13/
     +       30H OPTDRD    STEP TOLERANCE    =,E20.13/
     +       30H OPTDRD    GRADIENT TOLERANCE=,E20.13/
     +       30H OPTDRD    TRUST REG RADIUS  =,E20.13/
     +       30H OPTDRD    REL NOISE IN FCN  =,E20.13/
     +       30H OPTDRD    ANAL-FD TOLERANCE =,E20.13)
  906 FORMAT(52H OPTDRD    SHIFT FROM FORWARD TO CENTRAL DIFFERENCES,
     1   14H IN ITERATION , I5)
      END
      SUBROUTINE OPTF0D(NR,N,X,FCN,XPLS,FPLS,GPLS,ITRMCD,A,WRK)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C PROVIDE SIMPLEST INTERFACE TO MINIMIZATION PACKAGE.
C USER HAS NO CONTROL OVER OPTIONS.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> INITIAL ESTIMATE OF MINIMUM
C FCN          --> NAME OF ROUTINE TO EVALUATE MINIMIZATION FUNCTION.
C                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE.
C XPLS(N)     <--  LOCAL MINIMUM
C FPLS        <--  FUNCTION VALUE AT LOCAL MINIMUM XPLS
C GPLS(N)     <--  GRADIENT AT LOCAL MINIMUM XPLS
C ITRMCD      <--  TERMINATION CODE
C A(N,N)       --> WORKSPACE
C WRK(N,9)     --> WORKSPACE
C
C UNCMIN PACKAGE
C   CHANGES MADE:
C     NEW DRIVER:  UNCMND
C     PROBLEMS OF SIZE N=1 ALLOWED (SEE DFALTD "C---")
C
      DIMENSION X(N),XPLS(N),GPLS(N)
      DIMENSION A(NR,1),WRK(NR,1)
      EXTERNAL FCN,D1FCND,D2FCND
C
C EQUIVALENCE WRK(N,1) = UDIAG(N)
C             WRK(N,2) = G(N)
C             WRK(N,3) = P(N)
C             WRK(N,4) = TYPSIZ(N)
C             WRK(N,5) = SX(N)
C             WRK(N,6) = WRK0(N)
C             WRK(N,7) = WRK1(N)
C             WRK(N,8) = WRK2(N)
C             WRK(N,9) = WRK3(N)
C
      CALL DFALTD(N,X,WRK(1,4),FSCALE,METHOD,IEXP,MSG,NDIGIT,
     +     ITNLIM,IAGFLG,IAHFLG,IPR,DLT,GRADTL,STEPMX,STEPTL)
      CALL OPTDRD(NR,N,X,FCN,D1FCND,D2FCND,WRK(1,4),FSCALE,
     +     METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
     +     DLT,GRADTL,STEPMX,STEPTL,
     +     XPLS,FPLS,GPLS,ITRMCD,
     +     A,WRK(1,1),WRK(1,2),WRK(1,3),WRK(1,5),WRK(1,6),
     +     WRK(1,7),WRK(1,8),WRK(1,9))
      RETURN
      END
      SUBROUTINE OPTSTD(N,XPLS,FPLS,GPLS,X,ITNCNT,ICSCMX,
     +      ITRMCD,GRADTL,STEPTL,SX,FSCALE,ITNLIM,IRETCD,MXTAKE,IPR,MSG)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C UNCONSTRAINED MINIMIZATION STOPPING CRITERIA
C --------------------------------------------
C FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY
C OF THE FOLLOWING:
C 1) PROBLEM SOLVED WITHIN USER TOLERANCE
C 2) CONVERGENCE WITHIN USER TOLERANCE
C 3) ITERATION LIMIT REACHED
C 4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED
C
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE X[K]
C FPLS         --> FUNCTION VALUE AT NEW ITERATE F(XPLS)
C GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE
C X(N)         --> OLD ITERATE X[K-1]
C ITNCNT       --> CURRENT ITERATION K
C ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C ITRMCD      <--  TERMINATION CODE
C GRADTL       --> TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED CLOSE
C                  ENOUGH TO ZERO TO TERMINATE ALGORITHM
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION
C ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS
C IRETCD       --> RETURN CODE
C MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C MSG          --> IF MSG INCLUDES A TERM 8, SUPPRESS OUTPUT
C
C
      INTEGER N,ITNCNT,ICSCMX,ITRMCD,ITNLIM
      DIMENSION SX(N)
      DIMENSION XPLS(N),GPLS(N),X(N)
      LOGICAL MXTAKE
C
      ITRMCD=0
C
C LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X
      IF(IRETCD.NE.1) GO TO 50
C     IF(IRETCD.EQ.1)
C     THEN
        JTRMCD=3
        GO TO 600
C     ENDIF
   50 CONTINUE
C
C FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM.
C CHECK WHETHER WITHIN TOLERANCE
C
      D=MAX(ABS(FPLS),FSCALE)
      RGX=0.0D0
      DO 100 I=1,N
        RELGRD=ABS(GPLS(I))*MAX(ABS(XPLS(I)),1.D0/SX(I))/D
        RGX=MAX(RGX,RELGRD)
  100 CONTINUE
      JTRMCD=1
      IF(RGX.LE.GRADTL) GO TO 600
C
      IF(ITNCNT.EQ.0) RETURN
C
C FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM
C CHECK WHETHER WITHIN TOLERANCE.
C
      RSX=0.0D0
      DO 120 I=1,N
        RELSTP=ABS(XPLS(I)-X(I))/MAX(ABS(XPLS(I)),1.D0/SX(I))
        RSX=MAX(RSX,RELSTP)
  120 CONTINUE
      JTRMCD=2
      IF(RSX.LE.STEPTL) GO TO 600
C
C CHECK ITERATION LIMIT
C
      JTRMCD=4
      IF(ITNCNT.GE.ITNLIM) GO TO 600
C
C CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX
C
      IF(MXTAKE) GO TO 140
C     IF(.NOT.MXTAKE)
C     THEN
        ICSCMX=0
        RETURN
C     ELSE
  140   CONTINUE
        IF (MOD(MSG/8,2) .EQ. 0) WRITE(IPR,900)
        ICSCMX=ICSCMX+1
        IF(ICSCMX.LT.5) RETURN
        JTRMCD=5
C     ENDIF
C
C
C PRINT TERMINATION CODE
C
  600 ITRMCD=JTRMCD
      IF (MOD(MSG/8,2) .EQ. 0) GO TO(601,602,603,604,605), ITRMCD
      GO TO 700
  601 WRITE(IPR,901)
      GO TO 700
  602 WRITE(IPR,902)
      GO TO 700
  603 WRITE(IPR,903)
      GO TO 700
  604 WRITE(IPR,904)
      GO TO 700
  605 WRITE(IPR,905)
C
  700 RETURN
C
  900 FORMAT(48H0OPTSTD    STEP OF MAXIMUM LENGTH (STEPMX) TAKEN)
  901 FORMAT(43H0OPTSTD    RELATIVE GRADIENT CLOSE TO ZERO./
     +       48H OPTSTD    CURRENT ITERATE IS PROBABLY SOLUTION.)
  902 FORMAT(48H0OPTSTD    SUCCESSIVE ITERATES WITHIN TOLERANCE./
     +       48H OPTSTD    CURRENT ITERATE IS PROBABLY SOLUTION.)
  903 FORMAT(52H0OPTSTD    LAST GLOBAL STEP FAILED TO LOCATE A POINT,
     +       14H LOWER THAN X./
     +       51H OPTSTD    EITHER X IS AN APPROXIMATE LOCAL MINIMUM,
     +       17H OF THE FUNCTION,/
     +       50H OPTSTD    THE FUNCTION IS TOO NON-LINEAR FOR THIS,
     +       11H ALGORITHM,/
     +       34H OPTSTD    OR STEPTL IS TOO LARGE.)
  904 FORMAT(36H0OPTSTD    ITERATION LIMIT EXCEEDED./
     +       28H OPTSTD    ALGORITHM FAILED.)
  905 FORMAT(39H0OPTSTD    MAXIMUM STEP SIZE EXCEEDED 5,
     +       19H CONSECUTIVE TIMES./
     +       50H OPTSTD    EITHER THE FUNCTION IS UNBOUNDED BELOW,/
     +       47H OPTSTD    BECOMES ASYMPTOTIC TO A FINITE VALUE,
     +       30H FROM ABOVE IN SOME DIRECTION,/
     +       33H OPTSTD    OR STEPMX IS TOO SMALL)
      END
      SUBROUTINE QRAX1D(NR,N,R,I)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C INTERCHANGE ROWS I,I+1 OF THE UPPER HESSENBERG MATRIX R,
C COLUMNS I TO N
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF MATRIX
C R(N,N)      <--> UPPER HESSENBERG MATRIX
C I            --> INDEX OF ROW TO INTERCHANGE (I.LT.N)
C
      DIMENSION R(NR,1)
      DO 10 J=I,N
        TMP=R(I,J)
        R(I,J)=R(I+1,J)
        R(I+1,J)=TMP
   10 CONTINUE
      RETURN
      END
      SUBROUTINE QRAX2D(NR,N,R,I,A,B)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C PRE-MULTIPLY R BY THE JACOBI ROTATION J(I,I+1,A,B)
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF MATRIX
C R(N,N)      <--> UPPER HESSENBERG MATRIX
C I            --> INDEX OF ROW
C A            --> SCALAR
C B            --> SCALAR
C
      DIMENSION R(NR,1)
      DEN=SQRT(A*A + B*B)
      C=A/DEN
      S=B/DEN
      DO 10 J=I,N
        Y=R(I,J)
        Z=R(I+1,J)
        R(I,J)=C*Y - S*Z
        R(I+1,J)=S*Y + C*Z
   10 CONTINUE
      RETURN
      END
      SUBROUTINE QRUPDD(NR,N,A,U,V)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C FIND AN ORTHOGONAL (N*N) MATRIX (Q*) AND AN UPPER TRIANGULAR (N*N)
C MATRIX (R*) SUCH THAT (Q*)(R*)=R+U(V+)
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C A(N,N)      <--> ON INPUT:  CONTAINS R
C                  ON OUTPUT: CONTAINS (R*)
C U(N)         --> VECTOR
C V(N)         --> VECTOR
C
      DIMENSION A(NR,1)
      DIMENSION U(N),V(N)
C
C DETERMINE LAST NON-ZERO IN U(.)
C
      K=N
   10 IF(U(K).NE.0.D0 .OR. K.EQ.1) GO TO 20
C     IF(U(K).EQ.0.D0 .AND. K.GT.1)
C     THEN
        K=K-1
        GO TO 10
C     ENDIF
C
C (K-1) JACOBI ROTATIONS TRANSFORM
C     R + U(V+) --> (R*) + (U(1)*E1)(V+)
C WHICH IS UPPER HESSENBERG
C
   20 IF(K.LE.1) GO TO 40
        KM1=K-1
        DO 30 II=1,KM1
          I=KM1-II+1
          IF(U(I).NE.0.D0) GO TO 25
C         IF(U(I).EQ.0.)
C         THEN
            CALL QRAX1D(NR,N,A,I)
            U(I)=U(I+1)
            GO TO 30
C         ELSE
   25       CALL QRAX2D(NR,N,A,I,U(I),-U(I+1))
            U(I)=SQRT(U(I)*U(I) + U(I+1)*U(I+1))
C         ENDIF
   30   CONTINUE
C     ENDIF
C
C R <-- R + (U(1)*E1)(V+)
C
   40 DO 50 J=1,N
        A(1,J)=A(1,J) +U(1)*V(J)
   50 CONTINUE
C
C (K-1) JACOBI ROTATIONS TRANSFORM UPPER HESSENBERG R
C TO UPPER TRIANGULAR (R*)
C
      IF(K.LE.1) GO TO 100
        KM1=K-1
        DO 80 I=1,KM1
          IF(A(I,I).NE.0.D0) GO TO 70
C         IF(A(I,I).EQ.0.)
C         THEN
            CALL QRAX1D(NR,N,A,I)
            GO TO 80
C         ELSE
   70       T1=A(I,I)
            T2=-A(I+1,I)
            CALL QRAX2D(NR,N,A,I,T1,T2)
C         ENDIF
   80   CONTINUE
C     ENDIF
  100 RETURN
      END
      REAL FUNCTION R1MACH(I) 

c*********************************************************************72
c
C***BEGIN PROLOGUE  R1MACH
C***DATE WRITTEN   790101   (YYMMDD)
C***REVISION DATE  831014   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS) 
C***PURPOSE  Returns single precision machine dependent constants
C***DESCRIPTION
C     From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C
C     R1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function 
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          A = R1MACH(I)
C
C     where I=1,...,5.  The (output) value of A above is
C     determined by the (input) value of I.  The results for 
C     various values of I are discussed below.
C
C  Single-Precision Machine Constants
C  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. 
C  R1MACH(3) = B**(-T), the smallest relative spacing.
C  R1MACH(4) = B**(1-T), the largest relative spacing. 
C  R1MACH(5) = LOG10(B)
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5) 
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
C
C      DATA RMACH(1) / O"00014000000000000000" /
C      DATA RMACH(2) / O"37767777777777777777" /
C      DATA RMACH(3) / O"16404000000000000000" /
C      DATA RMACH(4) / O"16414000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
C
C     DATA RMACH(1) / X'9000400000000000' /
C     DATA RMACH(2) / X'6FFF7FFFFFFFFFFF' /
C     DATA RMACH(3) / X'FFA3400000000000' /
C     DATA RMACH(4) / X'FFA4400000000000' /
C     DATA RMACH(5) / X'FFD04D104D427DE8' /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
C
      DATA RMACH/1.18E-38,3.40E+38,0.595E-07,1.19E-07,0.30102999566/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C
C     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
C       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
C       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
C
C     DATA SMALL(1) / X'00800000' /
C     DATA LARGE(1) / X'7F7FFFFF' /
C     DATA RIGHT(1) / X'33800000' /
C     DATA DIVER(1) / X'34000000' /
C     DATA LOG10(1) / X'3E9A209B' /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C    (EXPRESSED IN INTEGER AND HEXADECIMAL)
C  *** THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS***
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C  ***THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS***
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1  .OR.  I .GT. 5)
     1   CALL XERROR ( 'R1MACH -- I OUT OF BOUNDS',25,1,2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END 
      SUBROUTINE RESLTD(NR,N,X,F,G,A,P,ITNCNT,IFLG,IPR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C PRINT INFORMATION
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> ITERATE X[K]
C F            --> FUNCTION VALUE AT X[K]
C G(N)         --> GRADIENT AT X[K]
C A(N,N)       --> HESSIAN AT X[K]
C P(N)         --> STEP TAKEN
C ITNCNT       --> ITERATION NUMBER K
C IFLG         --> FLAG CONTROLLING INFO TO PRINT
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C
      DIMENSION X(N),G(N),P(N),A(NR,1)
C PRINT ITERATION NUMBER
      WRITE(IPR,903) ITNCNT
      IF(IFLG.EQ.0) GO TO 120
C
C PRINT STEP
      WRITE(IPR,907)
      WRITE(IPR,905) (P(I),I=1,N)
C
C PRINT CURRENT ITERATE
  120 CONTINUE
      WRITE(IPR,904)
      WRITE(IPR,905) (X(I),I=1,N)
C
C PRINT FUNCTION VALUE
      WRITE(IPR,906)
      WRITE(IPR,905) F
C
C PRINT GRADIENT
      WRITE(IPR,908)
      WRITE(IPR,905) (G(I),I=1,N)
C
C PRINT HESSIAN FROM ITERATION K
      IF(IFLG.EQ.0) GO TO 140
      WRITE(IPR,901)
      DO 130 I=1,N
        WRITE(IPR,900) I
        WRITE(IPR,902) (A(I,J),J=1,I)
  130 CONTINUE
C
  140 RETURN
  900 FORMAT(15H RESLTD     ROW,I5)
  901 FORMAT(29H RESLTD       HESSIAN AT X(K))
  902 FORMAT(14H RESLTD       ,5(2X,E20.13))
  903 FORMAT(/21H0RESLTD    ITERATE K=,I5)
  904 FORMAT(18H RESLTD       X(K))
  905 FORMAT(22H RESLTD               ,5(2X,E20.13) )
  906 FORMAT(30H RESLTD       FUNCTION AT X(K))
  907 FORMAT(18H RESLTD       STEP)
  908 FORMAT(30H RESLTD       GRADIENT AT X(K))
      END
      SUBROUTINE SCLMLD(N,S,V,Z)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C MULTIPLY VECTOR BY SCALAR
C RESULT VECTOR MAY BE OPERAND VECTOR
C
C PARAMETERS
C ----------
C N            --> DIMENSION OF VECTORS
C S            --> SCALAR
C V(N)         --> OPERAND VECTOR
C Z(N)        <--  RESULT VECTOR
      DIMENSION V(N),Z(N)
      DO 100 I=1,N
        Z(I)=S*V(I)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SECFCD(NR,N,X,G,A,XPLS,GPLS,EPSM,ITNCNT,RNF,
     +     IAGFLG,NOUPDT,S,Y,U,W)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C UPDATE HESSIAN BY THE BFGS FACTORED METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE, X[K-1]
C G(N)         --> GRADIENT OR APPROXIMATE AT OLD ITERATE
C A(N,N)      <--> ON ENTRY: CHOLESKY DECOMPOSITION OF HESSIAN IN
C                    LOWER PART AND DIAGONAL.
C                  ON EXIT:  UPDATED CHOLESKY DECOMPOSITION OF HESSIAN
C                    IN LOWER TRIANGULAR PART AND DIAGONAL
C XPLS(N)      --> NEW ITERATE, X[K]
C GPLS(N)      --> GRADIENT OR APPROXIMATE AT NEW ITERATE
C EPSM         --> MACHINE EPSILON
C ITNCNT       --> ITERATION COUNT
C RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED, =0 ITHERWISE
C NOUPDT      <--> BOOLEAN: NO UPDATE YET
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C S(N)         --> WORKSPACE
C Y(N)         --> WORKSPACE
C U(N)         --> WORKSPACE
C W(N)         --> WORKSPACE
C
      DIMENSION X(N),XPLS(N),G(N),GPLS(N)
      DIMENSION A(NR,1)
      DIMENSION S(N),Y(N),U(N),W(N)
      LOGICAL NOUPDT,SKPUPD
C
      IF(ITNCNT.EQ.1) NOUPDT=.TRUE.
      DO 10 I=1,N
        S(I)=XPLS(I)-X(I)
        Y(I)=GPLS(I)-G(I)
   10 CONTINUE
      DEN1=DDOT(N,S,1,Y,1)
      SNORM2=DNRM2(N,S,1)
      YNRM2=DNRM2(N,Y,1)
      IF(DEN1.LT.SQRT(EPSM)*SNORM2*YNRM2) GO TO 110
C     IF(DEN1.GE.SQRT(EPSM)*SNORM2*YNRM2)
C     THEN
        CALL MVMLUD(NR,N,A,S,U)
        DEN2=DDOT(N,U,1,U,1)
C
C       L <-- SQRT(DEN1/DEN2)*L
C
        ALP=SQRT(DEN1/DEN2)
        IF(.NOT.NOUPDT) GO TO 50
C       IF(NOUPDT)
C       THEN
          DO 30 J=1,N
            U(J)=ALP*U(J)
            DO 20 I=J,N
              A(I,J)=ALP*A(I,J)
   20       CONTINUE
   30     CONTINUE
          NOUPDT=.FALSE.
          DEN2=DEN1
          ALP=1.0D0
C       ENDIF
   50   SKPUPD=.TRUE.
C
C       W = L(L+)S = HS
C
        CALL MVMLLD(NR,N,A,U,W)
        I=1
        IF(IAGFLG.NE.0) GO TO 55
C       IF(IAGFLG.EQ.0)
C       THEN
          RELTOL=SQRT(RNF)
          GO TO 60
C       ELSE
   55     RELTOL=RNF
C       ENDIF
   60   IF(I.GT.N .OR. .NOT.SKPUPD) GO TO 70
C       IF(I.LE.N .AND. SKPUPD)
C       THEN
          IF(ABS(Y(I)-W(I)) .LT. RELTOL*MAX(ABS(G(I)),ABS(GPLS(I))))
     +         GO TO 65
C         IF(ABS(Y(I)-W(I)) .GE. RELTOL*AMAX1(ABS(G(I)),ABS(GPLS(I))))
C         THEN
            SKPUPD=.FALSE.
            GO TO 60
C         ELSE
   65       I=I+1
            GO TO 60
C         ENDIF
C       ENDIF
   70   IF(SKPUPD) GO TO 110
C       IF(.NOT.SKPUPD)
C       THEN
C
C         W=Y-ALP*L(L+)S
C
          DO 75 I=1,N
            W(I)=Y(I)-ALP*W(I)
   75     CONTINUE
C
C         ALP=1/SQRT(DEN1*DEN2)
C
          ALP=ALP/DEN1
C
C         U=(L+)/SQRT(DEN1*DEN2) = (L+)S/SQRT((Y+)S * (S+)L(L+)S)
C
          DO 80 I=1,N
            U(I)=ALP*U(I)
   80     CONTINUE
C
C         COPY L INTO UPPER TRIANGULAR PART.  ZERO L.
C
          IF(N.EQ.1) GO TO 93
          DO 90 I=2,N
            IM1=I-1
            DO 85 J=1,IM1
              A(J,I)=A(I,J)
              A(I,J)=0.D0
   85       CONTINUE
   90     CONTINUE
C
C         FIND Q, (L+) SUCH THAT  Q(L+) = (L+) + U(W+)
C
   93     CALL QRUPDD(NR,N,A,U,W)
C
C         UPPER TRIANGULAR PART AND DIAGONAL OF A NOW CONTAIN UPDATED
C         CHOLESKY DECOMPOSITION OF HESSIAN.  COPY BACK TO LOWER
C         TRIANGULAR PART.
C
          IF(N.EQ.1) GO TO 110
          DO 100 I=2,N
            IM1=I-1
            DO 95 J=1,IM1
              A(I,J)=A(J,I)
   95       CONTINUE
  100     CONTINUE
C       ENDIF
C     ENDIF
  110 RETURN
      END
      SUBROUTINE SECNFD(NR,N,X,G,A,UDIAG,XPLS,GPLS,EPSM,ITNCNT,
     +     RNF,IAGFLG,NOUPDT,S,Y,T)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C UPDATE HESSIAN BY THE BFGS UNFACTORED METHOD
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE, X[K-1]
C G(N)         --> GRADIENT OR APPROXIMATE AT OLD ITERATE
C A(N,N)      <--> ON ENTRY: APPROXIMATE HESSIAN AT OLD ITERATE
C                    IN UPPER TRIANGULAR PART (AND UDIAG)
C                  ON EXIT:  UPDATED APPROX HESSIAN AT NEW ITERATE
C                    IN LOWER TRIANGULAR PART AND DIAGONAL
C                  [LOWER TRIANGULAR PART OF SYMMETRIC MATRIX]
C UDIAG        --> ON ENTRY: DIAGONAL OF HESSIAN
C XPLS(N)      --> NEW ITERATE, X[K]
C GPLS(N)      --> GRADIENT OR APPROXIMATE AT NEW ITERATE
C EPSM         --> MACHINE EPSILON
C ITNCNT       --> ITERATION COUNT
C RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN
C IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED, =0 OTHERWISE
C NOUPDT      <--> BOOLEAN: NO UPDATE YET
C                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C S(N)         --> WORKSPACE
C Y(N)         --> WORKSPACE
C T(N)         --> WORKSPACE
C
      DIMENSION X(N),G(N),XPLS(N),GPLS(N)
      DIMENSION A(NR,1)
      DIMENSION UDIAG(N)
      DIMENSION S(N),Y(N),T(N)
      LOGICAL NOUPDT,SKPUPD
C
C COPY HESSIAN IN UPPER TRIANGULAR PART AND UDIAG TO
C LOWER TRIANGULAR PART AND DIAGONAL
C
      DO 5 J=1,N
        A(J,J)=UDIAG(J)
        IF(J.EQ.N) GO TO 5
        JP1=J+1
        DO 4 I=JP1,N
          A(I,J)=A(J,I)
    4   CONTINUE
    5 CONTINUE
C
      IF(ITNCNT.EQ.1) NOUPDT=.TRUE.
      DO 10 I=1,N
        S(I)=XPLS(I)-X(I)
        Y(I)=GPLS(I)-G(I)
   10 CONTINUE
      DEN1=DDOT(N,S,1,Y,1)
      SNORM2=DNRM2(N,S,1)
      YNRM2=DNRM2(N,Y,1)
      IF(DEN1.LT.SQRT(EPSM)*SNORM2*YNRM2) GO TO 100
C     IF(DEN1.GE.SQRT(EPSM)*SNORM2*YNRM2)
C     THEN
        CALL MVMLSD(NR,N,A,S,T)
        DEN2=DDOT(N,S,1,T,1)
        IF(.NOT. NOUPDT) GO TO 50
C       IF(NOUPDT)
C       THEN
C
C         H <-- [(S+)Y/(S+)HS]H
C
          GAM=DEN1/DEN2
          DEN2=GAM*DEN2
          DO 30 J=1,N
            T(J)=GAM*T(J)
            DO 20 I=J,N
              A(I,J)=GAM*A(I,J)
   20       CONTINUE
   30     CONTINUE
          NOUPDT=.FALSE.
C       ENDIF
   50   SKPUPD=.TRUE.
C
C       CHECK UPDATE CONDITION ON ROW I
C
        DO 60 I=1,N
          TOL=RNF*MAX(ABS(G(I)),ABS(GPLS(I)))
          IF(IAGFLG.EQ.0) TOL=TOL/SQRT(RNF)
          IF(ABS(Y(I)-T(I)).LT.TOL) GO TO 60
C         IF(ABS(Y(I)-T(I)).GE.TOL)
C         THEN
            SKPUPD=.FALSE.
            GO TO 70
C         ENDIF
   60   CONTINUE
   70   IF(SKPUPD) GO TO 100
C       IF(.NOT.SKPUPD)
C       THEN
C
C         BFGS UPDATE
C
          DO 90 J=1,N
            DO 80 I=J,N
              A(I,J)=A(I,J)+Y(I)*Y(J)/DEN1-T(I)*T(J)/DEN2
   80       CONTINUE
   90     CONTINUE
C       ENDIF
C     ENDIF
  100 RETURN
      END
      SUBROUTINE SNDFDD(NR,N,XPLS,FCN,FPLS,A,SX,RNOISE,STEPSZ,ANBR)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C PURPOSE
C -------
C FIND SECOND ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A"
C TO THE SECOND DERIVATIVE (HESSIAN) OF THE FUNCTION DEFINED BY THE SUBP
C "FCN" EVALUATED AT THE NEW ITERATE "XPLS"
C
C FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE
C 1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION
C    IF NO ANALYTICAL USER FUNCTION HAS BEEN SUPPLIED FOR EITHER
C    THE GRADIENT OR THE HESSIAN AND IF THE OPTIMIZATION FUNCTION
C    "FCN" IS INEXPENSIVE TO EVALUATE.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C XPLS(N)      --> NEW ITERATE:   X[K]
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C FPLS         --> FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C A(N,N)      <--  FINITE DIFFERENCE APPROXIMATION TO HESSIAN
C                  ONLY LOWER TRIANGULAR MATRIX AND DIAGONAL
C                  ARE RETURNED
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C RNOISE       --> RELATIVE NOISE IN FNAME [F(X)]
C STEPSZ(N)    --> WORKSPACE (STEPSIZE IN I-TH COMPONENT DIRECTION)
C ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION)
C
C
      DIMENSION XPLS(N)
      DIMENSION SX(N)
      DIMENSION STEPSZ(N),ANBR(N)
      DIMENSION A(NR,1)
C
C FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION
C OF I-TH UNIT VECTOR.
C
      OV3 = 1.0D0/3.0D0
      DO 10 I=1,N
        STEPSZ(I)=RNOISE**OV3 * MAX(ABS(XPLS(I)),1.D0/SX(I))
        XTMPI=XPLS(I)
        XPLS(I)=XTMPI+STEPSZ(I)
        CALL FCN(N,XPLS,ANBR(I))
        XPLS(I)=XTMPI
   10 CONTINUE
C
C CALCULATE COLUMN I OF A
C
      DO 30 I=1,N
        XTMPI=XPLS(I)
        XPLS(I)=XTMPI+2.0D0*STEPSZ(I)
        CALL FCN(N,XPLS,FHAT)
        A(I,I)=((FPLS-ANBR(I))+(FHAT-ANBR(I)))/(STEPSZ(I)*STEPSZ(I))
C
C CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN
        IF(I.EQ.N) GO TO 25
        XPLS(I)=XTMPI+STEPSZ(I)
        IP1=I+1
        DO 20 J=IP1,N
          XTMPJ=XPLS(J)
          XPLS(J)=XTMPJ+STEPSZ(J)
          CALL FCN(N,XPLS,FHAT)
          A(J,I)=((FPLS-ANBR(I))+(FHAT-ANBR(J)))/(STEPSZ(I)*STEPSZ(J))
          XPLS(J)=XTMPJ
   20   CONTINUE
   25   XPLS(I)=XTMPI
   30 CONTINUE
      RETURN
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
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

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      SUBROUTINE TRGUPD(NR,N,X,F,G,A,FCN,SC,SX,NWTAKE,STEPMX,STEPTL,
     +     DLT,IRETCD,XPLSP,FPLSP,XPLS,FPLS,MXTAKE,IPR,METHOD,UDIAG)

c*********************************************************************72
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C PURPOSE
C -------
C DECIDE WHETHER TO ACCEPT XPLS=X+SC AS THE NEXT ITERATE AND UPDATE THE
C TRUST REGION DLT.
C
C PARAMETERS
C ----------
C NR           --> ROW DIMENSION OF MATRIX
C N            --> DIMENSION OF PROBLEM
C X(N)         --> OLD ITERATE X[K-1]
C F            --> FUNCTION VALUE AT OLD ITERATE, F(X)
C G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE
C A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN
C                  LOWER TRIANGULAR PART AND DIAGONAL.
C                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART
C FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION
C SC(N)        --> CURRENT STEP
C SX(N)        --> DIAGONAL SCALING MATRIX FOR X
C NWTAKE       --> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN
C STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE
C STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES
C                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM
C DLT         <--> TRUST REGION RADIUS
C IRETCD      <--> RETURN CODE
C                    =0 XPLS ACCEPTED AS NEXT ITERATE;
C                       DLT TRUST REGION FOR NEXT ITERATION.
C                    =1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE
C                       BECAUSE XPLS-X .LT. SMALLEST ALLOWABLE
C                       STEP LENGTH.
C                    =2 F(XPLS) TOO LARGE.  CONTINUE CURRENT ITERATION
C                       WITH NEW REDUCED DLT.
C                    =3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL
C                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO CONTINUE
C                       CURRENT ITERATION WITH NEW DOUBLED DLT.
C XPLSP(N)    <--> WORKSPACE [VALUE NEEDS TO BE RETAINED BETWEEN
C                  SUCCESIVE CALLS OF K-TH GLOBAL STEP]
C FPLSP       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]
C XPLS(N)     <--  NEW ITERATE X[K]
C FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)
C MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED
C IPR          --> DEVICE TO WHICH TO SEND OUTPUT
C METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM
C                    =1 LINE SEARCH
C                    =2 DOUBLE DOGLEG
C                    =3 MORE-HEBDON
C UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.)
C
      DIMENSION X(N),XPLS(N),G(N)
      DIMENSION SX(N),SC(N),XPLSP(N)
      DIMENSION A(NR,1)
      LOGICAL NWTAKE,MXTAKE
      DIMENSION UDIAG(N)
C
      IPR=IPR
      MXTAKE=.FALSE.
      DO 100 I=1,N
        XPLS(I)=X(I)+SC(I)
  100 CONTINUE
      CALL FCN(N,XPLS,FPLS)
      DLTF=FPLS-F
      SLP=DDOT(N,G,1,SC,1)
C
C NEXT STATEMENT ADDED FOR CASE OF COMPILERS WHICH DO NOT OPTIMIZE
C EVALUATION OF NEXT "IF" STATEMENT (IN WHICH CASE FPLSP COULD BE
C UNDEFINED).
      IF(IRETCD.EQ.4) FPLSP=0.0D0
C$    WRITE(IPR,961) IRETCD,FPLS,FPLSP,DLTF,SLP
      IF(IRETCD.NE.3 .OR. (FPLS.LT.FPLSP .AND. DLTF.LE. 1.D-4*SLP))
     +                                                     GO TO 130
C     IF(IRETCD.EQ.3 .AND. (FPLS.GE.FPLSP .OR. DLTF.GT. 1.E-4*SLP))
C     THEN
C
C       RESET XPLS TO XPLSP AND TERMINATE GLOBAL STEP
C
        IRETCD=0
        DO 110 I=1,N
          XPLS(I)=XPLSP(I)
  110   CONTINUE
        FPLS=FPLSP
        DLT=.5D0*DLT
C$      WRITE(IPR,951)
        GO TO 230
C     ELSE
C
C       FPLS TOO LARGE
C
  130   IF(DLTF.LE. 1.D-4*SLP) GO TO 170
C       IF(DLTF.GT. 1.E-4*SLP)
C       THEN
C$        WRITE(IPR,952)
          RLN=0.D0
          DO 140 I=1,N
            RLN=MAX(RLN,ABS(SC(I))/MAX(ABS(XPLS(I)),1.D0/SX(I)))
  140     CONTINUE
C$        WRITE(IPR,962) RLN
          IF(RLN.GE.STEPTL) GO TO 150
C         IF(RLN.LT.STEPTL)
C         THEN
C
C           CANNOT FIND SATISFACTORY XPLS SUFFICIENTLY DISTINCT FROM X
C
            IRETCD=1
C$          WRITE(IPR,954)
            GO TO 230
C         ELSE
C
C           REDUCE TRUST REGION AND CONTINUE GLOBAL STEP
C
  150       IRETCD=2
            DLTMP=-SLP*DLT/(2.D0*(DLTF-SLP))
C$          WRITE(IPR,963) DLTMP
            IF(DLTMP.GE. .1D0*DLT) GO TO 155
C           IF(DLTMP.LT. .1*DLT)
C           THEN
              DLT=.1D0*DLT
              GO TO 160
C           ELSE
  155         DLT=DLTMP
C           ENDIF
  160       CONTINUE
C$          WRITE(IPR,955)
            GO TO 230
C         ENDIF
C       ELSE
C
C         FPLS SUFFICIENTLY SMALL
C
  170     CONTINUE
C$        WRITE(IPR,958)
          DLTFP=0.D0
          IF (METHOD .EQ. 3) GO TO 180
C
C         IF (METHOD .EQ. 2)
C         THEN
C
          DO 177 I = 1, N
             TEMP = 0.0D0
             DO 173 J = I, N
                TEMP = TEMP + (A(J, I)*SC(J))
  173        CONTINUE
             DLTFP = DLTFP + TEMP*TEMP
  177     CONTINUE
          GO TO 190
C
C         ELSE
C
  180     DO 187 I = 1, N
             DLTFP = DLTFP + UDIAG(I)*SC(I)*SC(I)
             IF (I .EQ. N) GO TO 187
             TEMP = 0
             IP1 = I + 1
             DO 183 J = IP1, N
                TEMP = TEMP + A(I, J)*SC(I)*SC(J)
  183        CONTINUE
             DLTFP = DLTFP + 2.0D0*TEMP
  187     CONTINUE
C
C         END IF
C
  190     DLTFP = SLP + DLTFP/2.0D0
C$        WRITE(IPR,964) DLTFP,NWTAKE
          IF(IRETCD.EQ.2 .OR. (ABS(DLTFP-DLTF).GT. .1D0*ABS(DLTF))
     +         .OR. NWTAKE .OR. (DLT.GT. .99D0*STEPMX)) GO TO 210
C         IF(IRETCD.NE.2 .AND. (ABS(DLTFP-DLTF) .LE. .1*ABS(DLTF))
C    +         .AND. (.NOT.NWTAKE) .AND. (DLT.LE. .99*STEPMX))
C         THEN
C
C           DOUBLE TRUST REGION AND CONTINUE GLOBAL STEP
C
            IRETCD=3
            DO 200 I=1,N
              XPLSP(I)=XPLS(I)
  200       CONTINUE
            FPLSP=FPLS
            DLT=MIN(2.D0*DLT,STEPMX)
C$          WRITE(IPR,959)
            GO TO 230
C         ELSE
C
C           ACCEPT XPLS AS NEXT ITERATE.  CHOOSE NEW TRUST REGION.
C
  210       CONTINUE
C$          WRITE(IPR,960)
            IRETCD=0
            IF(DLT.GT. .99D0*STEPMX) MXTAKE=.TRUE.
            IF(DLTF.LT. .1D0*DLTFP) GO TO 220
C           IF(DLTF.GE. .1*DLTFP)
C           THEN
C
C             DECREASE TRUST REGION FOR NEXT ITERATION
C
              DLT=.5D0*DLT
              GO TO 230
C           ELSE
C
C             CHECK WHETHER TO INCREASE TRUST REGION FOR NEXT ITERATION
C
  220         IF(DLTF.LE. .75D0*DLTFP) DLT=MIN(2.D0*DLT,STEPMX)
C           ENDIF
C         ENDIF
C       ENDIF
C     ENDIF
  230 CONTINUE
C$    WRITE(IPR,953)
C$    WRITE(IPR,956) IRETCD,MXTAKE,DLT,FPLS
C$    WRITE(IPR,957)
C$    WRITE(IPR,965) (XPLS(I),I=1,N)
      RETURN
C
  951 FORMAT(55H TRGUPD    RESET XPLS TO XPLSP. TERMINATION GLOBAL STEP)
  952 FORMAT(26H TRGUPD    FPLS TOO LARGE.)
  953 FORMAT(38H0TRGUPD    VALUES AFTER CALL TO TRGUPD)
  954 FORMAT(54H TRGUPD    CANNOT FIND SATISFACTORY XPLS DISTINCT FROM,
     +       27H X.  TERMINATE GLOBAL STEP.)
  955 FORMAT(53H TRGUPD    REDUCE TRUST REGION. CONTINUE GLOBAL STEP.)
  956 FORMAT(21H TRGUPD       IRETCD=,I3/
     +       21H TRGUPD       MXTAKE=,L1/
     +       21H TRGUPD       DLT   =,E20.13/
     +       21H TRGUPD       FPLS  =,E20.13)
  957 FORMAT(32H TRGUPD       NEW ITERATE (XPLS))
  958 FORMAT(35H TRGUPD    FPLS SUFFICIENTLY SMALL.)
  959 FORMAT(54H TRGUPD    DOUBLE TRUST REGION.  CONTINUE GLOBAL STEP.)
  960 FORMAT(50H TRGUPD    ACCEPT XPLS AS NEW ITERATE.  CHOOSE NEW,
     +       38H TRUST REGION.  TERMINATE GLOBAL STEP.)
  961 FORMAT(18H TRGUPD    IRETCD=,I5/
     +       18H TRGUPD    FPLS  =,E20.13/
     +       18H TRGUPD    FPLSP =,E20.13/
     +       18H TRGUPD    DLTF  =,E20.13/
     +       18H TRGUPD    SLP   =,E20.13)
  962 FORMAT(18H TRGUPD    RLN   =,E20.13)
  963 FORMAT(18H TRGUPD    DLTMP =,E20.13)
  964 FORMAT(18H TRGUPD    DLTFP =,E20.13/
     +       18H TRGUPD    NWTAKE=,L1)
  965 FORMAT(14H TRGUPD       ,5(E20.13,3X))
      END
      SUBROUTINE UNCMND (N,X0,FCN,X,F,INFO,W,LW)

c*********************************************************************72
c
C***BEGIN PROLOGUE  UNCMND
C***DATE WRITTEN   870923    (YYMMDD)
C***REVISION DATE  871222    (YYMMDD)
C***CATEGORY NO.  G1B1A1
C***KEYWORDS  UNCONSTRAINED MINIMIZATION
C***AUTHOR  NASH, S.G., (GEORGE MASON UNIVERSITY)
C***PURPOSE  UNCMND minimizes a smooth nonlinear function of n variables.
C            A subroutine that computes the function value at any point
C            must be supplied, but derivative values are not required.
C            UNCMND provides a simple interface to more flexible lower
C            level routines.  User has no control over options.
C
C***DESCRIPTION
C     From the book, "Numerical Methods and Software" by
C                D. Kahaner, C. Moler, S. Nash
C                Prentice Hall, 1988
C
C     This routine uses a quasi-Newton algorithm with line search
C     to minimize the function represented by the subroutine FCN.
C     At each iteration, the nonlinear function is approximated
C     by a quadratic function derived from a Taylor series.
C     The quadratic function is minimized to obtain a search direction,
C     and an approximate minimum of the nonlinear function along
C     the search direction is found using a line search.  The
C     algorithm computes an approximation to the second derivative
C     matrix of the nonlinear function using quasi-Newton techniques.
C
C     The UNCMND package is quite general, and provides many options
C     for the user.  However, this subroutine is designed to be
C     easy to use, with few choices allowed.  For example:
C
C     1.  Only function values need be computed.  First derivative
C     values are obtained by finite-differencing.  This can be
C     very costly when the number of variables is large.
C
C     2.  It is assumed that the function values can be obtained
C     accurately (to an accuracy comparable to the precision of
C     the computer arithmetic).
C
C     3.  At most 150 iterations are allowed.
C
C     4.  It is assumed that the function values are well-scaled,
C     that is, that the optimal function value is not pathologically
C     large or small.
C
C     For more information, see the reference listed below.
C
C PARAMETERS
C ----------
C N            --> INTEGER
C                  Dimension of problem
C X0(N)        --> DOUBLE PRECISION
C                  Initial estimate of minimum
C FCN          --> Name of routine to evaluate minimization function.
C                  Must be declared EXTERNAL in calling routine, and
C                  have calling sequence
C                      SUBROUTINE FCN(N, X, F)
C                  with N and X as here, F the computed function value.
C X(N)        <--  DOUBLE PRECISION
C                  Local minimum
C F           <--  DOUBLE PRECISION
C                  Function value at local minimum X
C INFO        <--  INTEGER
C                  Termination code
C                      INFO =  0:  Optimal solution found
C                      INFO =  1:  Terminated with gradient small,
C                                  X is probably optimal
C                      INFO =  2:  Terminated with stepsize small,
C                                  X is probably optimal
C                      INFO =  3:  Lower point cannot be found,
C                                  X is probably optimal
C                      INFO =  4:  Iteration limit (150) exceeded
C                      INFO =  5:  Too many large steps,
C                                  function may be unbounded
C                      INFO = -1:  Insufficient workspace
C W(LW)        --> DOUBLE PRECISION
C                  Workspace
C LW           --> INTEGER
C                  Size of workspace, at least N*(N+10)
C
C***REFERENCES  R.B. SCHNABEL, J.E. KOONTZ, AND BE.E. WEISS, A MODULAR
C                 SYSTEM OF ALGORITHMS FOR UNCONSTRAINED MINIMIZATION,
C                 REPORT CU-CS-240-82, COMP. SCI. DEPT., UNIV. OF
C                 COLORADO AT BOULDER, 1982.
C***ROUTINES CALLED  OPTDRD, XERROR
C***END PROLOGUE  UNCMND
      IMPLICIT  DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(N),X(N),W(LW)
      CHARACTER ERRMSG*68
      EXTERNAL  FCN, D1FCND, D2FCND
c
C SUBDIVIDE WORKSPACE
c
C***FIRST EXECUTABLE STATEMENT  UNCMND
      IG  = 1
      IT  = IG  + N
      IW1 = IT  + N
      IW2 = IW1 + N
      IW3 = IW2 + N
      IW4 = IW3 + N
      IW5 = IW4 + N
      IW6 = IW5 + N
      IW7 = IW6 + N
      IW8 = IW7 + N
      IA  = IW8 + N
      LWMIN = IA + N*N-1
      IF (LWMIN .GT. LW) THEN
          INFO = -1
          WRITE(ERRMSG, '(
     *      ''UNCMND ERROR (INFO=-1) -- INSUFFICIENT WORKSPACE'',
     *      '', LW = '', I5 )' ) LW
          CALL XERROR(ERRMSG(1:60), 60, -1, 0)
          RETURN
      END IF
c
C SET UP PARAMETERS FOR OPTDRD
c
C PARAMETERS THAT SHOULD NOT BE CHANGED WHEN USING CONDENSED CODE
C
C NR     = PARAMETER USED TO DIVIDE WORKSPACE
C METHOD = 1 (LINE SEARCH) -- DO NOT CHANGE
C MSG    = 9 => NO PRINTING, N=1 ALLOWED
C IAGFLG = 1 => ANALYTIC GRADIENT SUPPLIED (0 OTHERWISE)
C IAHFLG = 1 => ANALYTIC HESSIAN  SUPPLIED (0 OTHERWISE)
C IPR    = DEVICE FOR OUTPUT (IRRELEVANT IN CURRENT VERSION)
C DLT    = (IRRELEVANT PARAMETER FOR METHOD = 1)
C EPSM   = MACHINE EPSILON
C IEXP   = 1 => FUNCTION EXPENSIVE TO EVALUATE (IEXP = 0 => CHEAP)
C
      NR = N
      METHOD = 1
      MSG = 9
      IAGFLG = 0
      IAHFLG = 0
      IPR = 0
      DLT = -1.0D0
      EPSM = D1MACH(4)
      IEXP = 1
C
C PARAMETERS THAT MAY BE CHANGED:
C
C NDIGIT = -1 => OPTDRD ASSUMES F IS FULLY ACCURATE
C ITNLIM = 150 = MAXIMUM NUMBER OF ITERATIONS ALLOWED
C GRADTL = ZERO TOLERANCE FOR GRADIENT, FOR CONVERGENCE TESTS
C STEPMX = MAXIMUM ALLOWABLE STEP SIZE
C STEPTL = ZERO TOLERANCE FOR STEP, FOR CONVERGENCE TESTS
C FSCALE = TYPICAL ORDER-OF-MAGNITUDE SIZE OF FUNCTION
C TYPSIZ = TYPICAL ORDER-OF-MAGNITUDE SIZE OF X (STORED IN W(LT))
C
      NDIGIT = -1
      ITNLIM = 150
      GRADTL = EPSM**(1.0D0/3.0D0)
      STEPMX = 0.0D0
      STEPTL = SQRT(EPSM)
      FSCALE = 1.0D0
      DO 10 LT = IT,IT+N-1
          W(LT) = 1.0D0
   10 CONTINUE
C
C MINIMIZE FUNCTION
C
      CALL OPTDRD (NR, N, X0, FCN, D1FCND, D2FCND, W(IT), FSCALE,
     +             METHOD, IEXP, MSG, NDIGIT, ITNLIM, IAGFLG, IAHFLG,
     +             IPR, DLT, GRADTL, STEPMX, STEPTL,
     +             X, F, W(IG), INFO, W(IA),
     +             W(IW1), W(IW2), W(IW3), W(IW4),
     +             W(IW5), W(IW6), W(IW7), W(IW8))
C
      IF (INFO .EQ. 1) THEN
          WRITE(ERRMSG, '(
     *      ''UNCMND WARNING -- INFO = 1'',
     *      '': PROBABLY CONVERGED, GRADIENT SMALL'')' )
          CALL XERROR(ERRMSG(1:62), 62, INFO, 0)
      END IF
      IF (INFO .EQ. 2) THEN
          WRITE(ERRMSG, '(
     *      ''UNCMND WARNING -- INFO = 2'',
     *      '': PROBABLY CONVERGED, STEPSIZE SMALL'')' )
          CALL XERROR(ERRMSG(1:62), 62, INFO, 0)
      END IF
      IF (INFO .EQ. 3) THEN
          WRITE(ERRMSG, '(
     *      ''UNCMND WARNING -- INFO = 3'',
     *      '': CANNOT FIND LOWER POINT'')' )
          CALL XERROR(ERRMSG(1:51), 51, INFO, 0)
      END IF
      IF (INFO .EQ. 4) THEN
          WRITE(ERRMSG, '(
     *      ''UNCMND WARNING -- INFO = 4'',
     *      '': TOO MANY ITERATIONS'')' )
          CALL XERROR(ERRMSG(1:47), 47, INFO, 0)
      END IF
      IF (INFO .EQ. 5) THEN
          WRITE(ERRMSG, '(
     *      ''UNCMND WARNING -- INFO = 5'',
     *      '': TOO MANY LARGE STEPS, POSSIBLY UNBOUNDED'')' )
          CALL XERROR(ERRMSG(1:68), 68, INFO, 0)
      END IF
C
      RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCLR

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERCLR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Resets current error number to zero.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        This routine simply resets the current error number to zero.
C        This may be necessary to do in order to determine that
C        a certain error has occurred again since the last time
C        NUMXER was referenced.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XERCLR
C***FIRST EXECUTABLE STATEMENT  XERCLR
      JUNK = J4SAVE(1,0,.TRUE.)
      RETURN
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERDMP

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERDMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints the error tables and then clears them.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERDMP prints the error tables, then clears them.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERSAV
C***END PROLOGUE  XERDMP
C***FIRST EXECUTABLE STATEMENT  XERDMP
      CALL XERSAV(' ',0,0,0,KOUNT)
      RETURN
      END
      SUBROUTINE XERMAX(MAX)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERMAX
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Sets maximum number of times any error message is to be
C            printed.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERMAX sets the maximum number of times any message
C        is to be printed.  That is, non-fatal messages are
C        not to be printed after they have occured MAX times.
C        Such non-fatal messages may be printed less than
C        MAX times even if they occur MAX times, if error
C        suppression mode (KONTRL=0) is ever in effect.
C
C     Description of Parameter
C      --Input--
C        MAX - the maximum number of times any one message
C              is to be printed.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XERMAX
C***FIRST EXECUTABLE STATEMENT  XERMAX
      JUNK = J4SAVE(4,MAX,.TRUE.)
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  870916   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  16 SEPT 1987
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
C         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            IF(IUNIT.EQ.0)THEN
              WRITE (*,'(1X,A)') MESSG(ICHAR:LAST)
            ELSE
              WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  870916   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  16 SEPT 1987
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
C            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) I1
                 IF (I.EQ.2) WRITE (*,FORM) I2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) I1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) I2
               ENDIF
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) R1
                 IF (I.EQ.2) WRITE (*,FORM) R2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) R1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) R2
               ENDIF
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               IF(IUINT.EQ.0)THEN
                 WRITE(*,30) LERR
               ELSE
                 WRITE (IUNIT,30) LERR
               ENDIF
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETF(KONTRL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XGETF
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns current value of error control flag.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C   Abstract
C        XGETF returns the current value of the error control flag
C        in KONTRL.  See subroutine XSETF for flag value meanings.
C        (KONTRL is an output parameter only.)
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETF
C***FIRST EXECUTABLE STATEMENT  XGETF
      KONTRL = J4SAVE(2,0,.FALSE.)
      RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE XGETUN(IUNIT)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XGETUN
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns the (first) output file to which messages are being
C            sent.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XGETUN gets the (first) output file to which error messages
C        are being sent.  To find out if more than one file is being
C        used, one must use the XGETUA routine.
C
C     Description of Parameter
C      --Output--
C        IUNIT - the logical unit number of the  (first) unit to
C                which error messages are being sent.
C                A value of zero means that the default file, as
C                defined by the I1MACH routine, is being used.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision --- 23 May 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUN
C***FIRST EXECUTABLE STATEMENT  XGETUN
      IUNIT = J4SAVE(3,0,.FALSE.)
      RETURN
      END
      SUBROUTINE XSETF(KONTRL)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XSETF
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3A
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Sets the error control flag.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XSETF sets the error control flag value to KONTRL.
C        (KONTRL is an input parameter only.)
C        The following table shows how each message is treated,
C        depending on the values of KONTRL and LEVEL.  (See XERROR
C        for description of LEVEL.)
C
C        If KONTRL is zero or negative, no information other than the
C        message itself (including numeric values, if any) will be
C        printed.  If KONTRL is positive, introductory messages,
C        trace-backs, etc., will be printed in addition to the message.
C
C              IABS(KONTRL)
C        LEVEL        0              1              2
C        value
C          2        fatal          fatal          fatal
C
C          1     not printed      printed         fatal
C
C          0     not printed      printed        printed
C
C         -1     not printed      printed        printed
C                                  only           only
C                                  once           once
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE,XERRWV
C***END PROLOGUE  XSETF
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF ((KONTRL.GE.(-2)).AND.(KONTRL.LE.2)) GO TO 10
         CALL XERRWV('XSETF  -- INVALID VALUE OF KONTRL (I1).',33,1,2,
     1  1,KONTRL,0,0,0.,0.)
         RETURN
   10 JUNK = J4SAVE(2,KONTRL,.TRUE.)
      RETURN
      END
      SUBROUTINE XSETUA(IUNITA,N)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XSETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3B
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Sets up to 5 unit numbers to which messages are to be sent.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XSETUA may be called to declare a list of up to five
C        logical units, each of which is to receive a copy of
C        each error message processed by this package.
C        The purpose of XSETUA is to allow simultaneous printing
C        of each error message on, say, a main output file,
C        an interactive terminal, and other files such as graphics
C        communication files.
C
C     Description of Parameters
C      --Input--
C        IUNIT - an array of up to five unit numbers.
C                Normally these numbers should all be different
C                (but duplicates are not prohibited.)
C        N     - the number of unit numbers provided in IUNIT
C                must have 1 .LE. N .LE. 5.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE,XERRWV
C***END PROLOGUE  XSETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XSETUA
      IF ((N.GE.1).AND.(N.LE.5)) GO TO 10
         CALL XERRWV('XSETUA -- INVALID VALUE OF N (I1).',34,1,2,
     1  1,N,0,0,0.,0.)
         RETURN
   10 CONTINUE
      DO 20 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         JUNK = J4SAVE(INDEX,IUNITA(I),.TRUE.)
   20 CONTINUE
      JUNK = J4SAVE(5,N,.TRUE.)
      RETURN
      END
      SUBROUTINE XSETUN(IUNIT)

c*********************************************************************72
c
C***BEGIN PROLOGUE  XSETUN
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3B
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Sets output file to which error messages are to be sent.
C***DESCRIPTION
C    From the book "Numericl Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XSETUN sets the output file to which error messages are to
C        be sent.  Only one file will be used.  See XSETUA for
C        how to declare more than one file.
C
C     Description of Parameter
C      --Input--
C        IUNIT - an input parameter giving the logical unit number
C                to which error messages are to be sent.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  7 June 1978
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XSETUN
C***FIRST EXECUTABLE STATEMENT  XSETUN
      JUNK = J4SAVE(3,IUNIT,.TRUE.)
      JUNK = J4SAVE(5,1,.TRUE.)
      RETURN
      END
