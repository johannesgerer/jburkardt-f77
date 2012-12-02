      SUBROUTINE DERPAR(N, X, XLOW, XUPP, EPS, W, INITAL, ITIN, HH,
     * HMAX, PREF, NDIR, E, MXADMS, NCORR, NCRAD, NOUT, OUT,
     * MAXOUT, NPRNT)

c*********************************************************************72
c
cc DERPAR carries out Kubicek's continuation algorithm.
c
      DIMENSION X(11), XLOW(11), XUPP(11), W(11), HMAX(11),
     * PREF(11), NDIR(11), OUT(100,12), F(11), G(10,11), BETA(11),
     * MARK(11), DXDT(11)
C OBTAINING OF DEPENDENCE OF SOLUTION X(ALFA) OF EQUATION
C F ( X , ALFA ) = 0
C ON PARAMETER ALFA BY MODIFIED METHOD OF DIFFERENTIATION
C WITH RESPECT TO PARAMETER
C N - NUMBER OF UNKNOWNS X(I)
C X(1),...,X(N) - INITIAL VALUES OF X(I), AFTER RETURN FINAL VALUES
C                 OF X(I)
C X(N+1) - INITIAL VALUE OF PARAMETER ALFA, AFTER RETURN FINAL VALUE
C          OF ALFA
C XLOW(1),...,XLOW(N) - LOWER BOUNDS FOR X(I)
C XUPP(1),...,XUPP(N) - UPPER BOUNDS FOR X(I)
C XLOW(N+1),XUPP(N+1) - LOWER AND UPPER BOUNDS FOR ALFA
C    IF XLOW OR XUPP IS EXCEEDED, THEN END OF DERPAR
C    AND MAXOUT=-2 AFTER RETURN
C EPS - ACCURACY DESIRED IN NEWTON ITERATION FOR
C    SUM OF (W(I)*ABS(XNEW(I)-XOLD(I))), I=1,...,N+1
C W(1),...,W(N+1) - WEIGHTS USED IN TERMINATION CRITERION OF NEWTON
C                   PROCESS
C INITAL - IF (INITAL.NE.0) THEN SEVERAL STEPS IN NEWTON ITERATION
C         ARE MADE BEFORE COMPUTATION IN ORDER TO INCREASE
C         ACCURACY OF INITIAL POINT-
C     IF (INITAL.EQ.1.AND.EPS-ACCURACY IS NOT FULFILLED IN ITIN
C     ITERATIONS) THEN RETURN.   IF (INITAL.EQ.2) THEN ALWAYS RETURN
C     AFTER INITIAL NEWTON ITERATION, RESULTS ARE IN X.
C     IF (INITAL.EQ.3) THEN CONTINUE IN DERPAR AFTER INITIAL NEWTON.
C     IF (INITAL.EQ.0) THEN NO INITIAL NEWTON ITERATION IS USED.
C ITIN - MAXIMAL NUMBER OF INITIAL NEWTON ITERATIONS. IF
C       EPS-ACCURACY IS NOT FULFILLED IN ITIN ITERATIONS THEN
C       ITIN=-1 AFTER RETURN.
C HH - INTEGRATION STEP ALONG ARC LENGTH OF SOLUTION LOCUS
C HMAX(1),...,HMAX(N+1) - UPPER BOUNDS FOR INCREMENTS OF X(I) IN
C                        ONE INTEGRATION STEP (APPROXIMATION ONLY)
C PREF(1),...,PREF(N+1) - PREFERENCE NUMBERS (EXPLANATION SEE IN
C                         SUBR.GAUSE)
C NDIR(1),...,NDIR(N+1) - INITIAL CHANGE OF X(I) IS POSITIVE ALONG
C     SOLUTION LOCUS (CURVE) IF NDIR(I)=1 AND NEGATIVE IF
C     NDIR(I)=-1.
C E - CRITERION FOR TEST ON CLOSED CURVE, IF
C    (SUM OF (W(I)*ABS(X(I)-XINITIAL(I))),I=1,...,N+1).LE.E)
C    THEN CLOSED CURVE MAY BE EXPECTED.
C MXADMS - MAXIMAL ORDER OF ADAMS-BASHFORTH FORMULA,
C          1.LE.MXADMS.LE.4.
C NCORR - MAXIMAL NUMBER OF NEWTON CORRECTIONS AFTER PREDICTION
C    BY ADAMS-BASHFORTH METHOD.
C NCRAD - IF (NCRAD.NE.0) THEN ADDITIONAL NEWTON CORRECTION
C    WITHOUT NEW COMPUTATION OF JACOBIAN MATRIX IS USED.
C NOUT - AFTER RETURN NUMBER OF CALCULATED POINTS ON THE CURVE
C        X(ALFA), NOUT.LE.MAXOUT.
C OUT(J,1),OUT(J,2),...,OUT(J,N+1) - J-TH POINT X(1),...,X(N),ALFA
C    ON CURVE X(ALFA)
C OUT(J,N+2) - VALUE OF SQRT(SUM OF SQUARES OF F).
C    IF (OUT(J,N+2).LT.0.0) THEN  ABS(OUT(J,N+2)) CORRESPONDS TO X
C    AND ALFA NOT EXACTLY, BECAUSE ADDITIONAL NEWTON CORRECTION WAS
C    USED (NCRAD.NE.0).
C    VALUES F(I) ARE NOT COMPUTED FOR X AND ALFA PRINTED AND/OR
C    STORED AND THEREFORE LAST TIME COMPUTED VALUE OF SQRT(SUM OF
C    SQUARES OF F) IS ON OUR DISPOSAL ONLY.
C MAXOUT - MAXIMAL NUMBER OF CALCULATED POINTS ON CURVE X(ALFA).
C    IF MAXOUT AFTER RETURN EQUALS TO-
C          -1 - THEN CLOSED CURVE X(ALFA) MAY BE EXPECTED
C          -2 - THEN BOUND XLOW OR XUPP WAS EXCEEDED
C          -3 - THEN SINGULAR JACOBIAN MATRIX OCCURRED, ITS RANK IS
C               .LT. N.
C NPRNT - IF (NPRNT.EQ.3) THEN RESULTING POINTS ON CURVE X(ALFA)
C    ARE IN ARRAY OUT( , ) AFTER RETURN. IF (NPRNT.EQ.1) THEN THESE
C    POINTS ARE PRINTED ONLY.  IF (NPRNT.EQ.2) THEN THESE POINTS ARE
C    BOTH PRINTED AND IN ARRAY OUT.
C SUBROUTINE FCTN MUST BE PROGRAMMED IN FOLLOWING FORM -
C   SUBROUTINE FCTN(N,X,F,G)
C   DIMENSION X(11),F(10),G(10,11)
C   F(I)= FI (X(1),X(2),...,X(N),ALFA) FOR I=1,...,N
C   G(I,J)= D FI (X(1),...,X(N),ALFA)/ D X(J) FOR  I,J=1,...,N
C   G(I,N+1)= D FI (X(1),...,X(N),ALFA)/ D ALFA  FOR I=1,...,N
C   RETURN
C   END
      DATA INDIC, INDSP /1H*,1H /
c
C LW IS PRINTER DEVICE NUMBER
c
      N1 = N + 1
      LW = 6

      IF (INITAL) 10, 60, 10
c
C INITIAL NEWTON ITERATIONS
c
   10 DO 40 L=1,ITIN
        CALL FCTN(N, X, F, G)
        SQUAR = 0.0
        DO 20 I=1,N
          SQUAR = SQUAR + F(I)**2
   20   CONTINUE
        LL = L - 1
        SQUAR = SQRT(SQUAR)
        IF (NPRNT.NE.3) WRITE (LW,99999) LL, (X(I),I=1,N1), SQUAR
        CALL GAUSE(N, G, F, M, 10, 11, PREF, BETA, K)
        IF (M.EQ.0) GO TO 310
        P = 0.0
        DO 30 J=1,N1
          X(J) = X(J) - F(J)
          P = P + ABS(F(J))*W(J)
   30   CONTINUE
        IF (P.LE.EPS) GO TO 50
   40 CONTINUE
      WRITE (LW,99998) ITIN
      ITIN = -1
      IF (INITAL.EQ.1) RETURN
   50 IF (NPRNT.NE.3) WRITE (LW,99997) (X(I),I=1,N1)
      IF (INITAL.EQ.2) RETURN
c
C AFTER INITIAL NEWTON ITERATIONS
c
   60 IF (NPRNT.NE.3) WRITE (LW,99996)
      KOUT = 0
      NOUT = 0
      MADMS = 0
      NC = 1
      K1 = 0
   70 CALL FCTN(N, X, F, G)
      SQUAR = 0.0
      DO 80 I=1,N
        SQUAR = SQUAR + F(I)**2
   80 CONTINUE
      CALL GAUSE(N, G, F, M, 10, 11, PREF, BETA, K)
      IF (M.EQ.0) GO TO 310
      IF (K1.EQ.K) GO TO 90
c
C CHANGE OF INDEPENDENT VARIABLE (ITS INDEX = K NOW)
c
      MADMS = 0
      K1 = K
   90 SQUAR = SQRT(SQUAR)
      IF (NCRAD.EQ.1) SQUAR = -SQUAR
      P = 0.0
      DO 100 I=1,N1
        P = P + W(I)*ABS(F(I))
  100 CONTINUE
      IF (P.LE.EPS) GO TO 130
      IF (NC.GE.NCORR) GO TO 120
      DO 110 I=1,N1
        X(I) = X(I) - F(I)
  110 CONTINUE
c
C ONE ITERATION IN NEWTON METHOD
c
      NC = NC + 1
      GO TO 70
  120 IF (NCORR.EQ.0) GO TO 130
      WRITE (LW,99995) NCORR, P
  130 NC = 1
      IF (NCRAD.EQ.0) GO TO 150
c
C ADDITIONAL NEWTON CORRECTION
c
      DO 140 I=1,N1
        X(I) = X(I) - F(I)
  140 CONTINUE
  150 NOUT = NOUT + 1
      DO 160 I=1,N1
        MARK(I) = INDSP
  160 CONTINUE
      MARK(K) = INDIC
      IF (NPRNT.EQ.3) GO TO 170
      WRITE (LW,99994)
      WRITE (LW,99993) (X(I),MARK(I),I=1,N1), SQUAR
  170 IF (NPRNT.EQ.1) GO TO 200
      IF (NOUT.LE.100) GO TO 180
      WRITE (LW,99992)
      RETURN
  180 DO 190 I=1,N1
        OUT(NOUT,I) = X(I)
  190 CONTINUE
      OUT(NOUT,N+2) = SQUAR
      GO TO 210
  200 IF (NOUT.EQ.1) GO TO 180
  210 IF (NOUT.GE.MAXOUT) RETURN
      DO 220 I=1,N1
        IF (X(I).LT.XLOW(I) .OR. X(I).GT.XUPP(I)) GO TO 300
  220 CONTINUE
      IF (NOUT.LE.3) GO TO 240
      P = 0.0
      DO 230 I=1,N1
        P = P + W(I)*ABS(X(I)-OUT(1,I))
  230 CONTINUE
      IF (P.LE.E) GO TO 290
c
C CLOSED CURVE MAY BE EXPECTED
c
  240 DXK2 = 1.0
      DO 250 I=1,N1
        DXK2 = DXK2 + BETA(I)**2
  250 CONTINUE
      DXDT(K) = 1.0/SQRT(DXK2)*FLOAT(NDIR(K))
c
C DERIVATIVE OF INDEPENDENT VARIABLE X(K) WITH RESPECT TO ARC LENGTH
C OF SOLUTION LOCUS IS COMPUTED HERE
c
      H = HH
      DO 270 I=1,N1
        NDIR(I) = 1
        IF (I.EQ.K) GO TO 260
        DXDT(I) = BETA(I)*DXDT(K)
  260   IF (DXDT(I).LT.0.0) NDIR(I) = -1
        IF (H*ABS(DXDT(I)).LE.HMAX(I)) GO TO 270
        MADMS = 0
        H = HMAX(I)/ABS(DXDT(I))
  270 CONTINUE
      IF (NOUT.LE.KOUT+3) GO TO 280
      IF (H*ABS(DXDT(K)).LE.0.8*ABS(X(K)-OUT(1,K))) GO TO 280
      IF ((OUT(1,K)-X(K))*FLOAT(NDIR(K)).LE.0.0) GO TO 280
      MADMS = 0
      IF (H*ABS(DXDT(K)).LE.ABS(X(K)-OUT(1,K))) GO TO 280
      H = ABS(X(K)-OUT(1,K))/ABS(DXDT(K))
      KOUT = NOUT
  280 CALL ADAMS(N, DXDT, MADMS, H, X, MXADMS)
      GO TO 70
  290 WRITE (LW,99991)
      MAXOUT = -1
      RETURN
  300 MAXOUT = -2
      RETURN
  310 WRITE (LW,99990) (X(I),I=1,N1)
      MAXOUT = -3
      RETURN
99999 FORMAT (3X, 7H DERPAR, I3, 25H.INITIAL NEWTON ITERATION/22X,
     * 11HX,ALFA,SQF=, 5F15.7/(33X, 5F15.7))
99998 FORMAT (/16H DERPAR OVERFLOW, I5, 2X, 18HINITIAL ITERATIONS/)
99997 FORMAT (3X, 39H DERPAR AFTER INITIAL NEWTON ITERATIONS/22X,
     * 7HX,ALFA=, 4X, 5F15.7/(33X, 5F15.7))
99996 FORMAT (/6X, 45H DERPAR RESULTS (VARIABLE CHOSEN AS INDEPENDE,
     * 18HNT IS MARKED BY *)/)
99995 FORMAT (/37H DERPAR NUMBER OF NEWTON CORRECTIONS=, I5,
     * 31H  IS NOT SUFFICIENT,ERROR OF X=, F15.7/)
99994 FORMAT (6X, 27H DERPAR RESULTS X,ALFA,SQF=)
99993 FORMAT (33X, 5(F15.7, A1))
99992 FORMAT (/28H DERPAR OUTPUT ARRAY IS FULL//)
99991 FORMAT (/36H DERPAR CLOSED CURVE MAY BE EXPECTED/)
99990 FORMAT (/48H DERPAR SINGULAR JACOBIAN MATRIX FOR X AND ALFA=,
     * 5F12.6/(48X, 5F12.6))
      end
      subroutine adams ( n, d, madms, h, x, mxadms )

c*********************************************************************72
c
cc ADAMS applies Adams-Bashforth methods to an ODE.
c
      real d(11)
      real der(4,11)
      real x(11)

      n1 = n + 1

      do i = 1, 3
        do j = 1, n + 1
          der(i+1,j) = der(i,j)
        end do
      end do

      madms = madms + 1
      if (madms.gt.mxadms) madms = mxadms
      if (madms.gt.4) madms = 4

      do i = 1, n + 1

        der(1,i) = d(i)

        if ( madms .eq. 1 ) then
          x(i) = x(i) + h*der(1,i)
        else if ( madms .eq. 2 ) then
          x(i) = x(i) + 0.5*h*(3.0*der(1,i)-der(2,i))
        else if ( madms .eq. 3 ) then
          x(i) = x(i) + h*(23.0*der(1,i)-16.0*der(2,i)+5.0*der(3,i))/
     &   12.0
        else if ( madms .eq. 4 ) then
          x(i) = x(i) + h*(55.0*der(1,i)-59.0*der(2,i)+37.0*der(3,i)
     &   -9.0*der(4,i))/24.0
        end if

      end do

      return
      end
      subroutine gause(n, a, b, m, nn, mm, pref, beta, k)

c*********************************************************************72
c
cc GAUSE solves a linear system of N equations in N+1 unknowns.
c
C SOLUTION OF N LINEAR EQUATIONS FOR N+1 UNKNOWNS
C BASED ON GAUSSIAN ELIMINATION WITH PIVOTING.
C N - NUMBER OF EQUATIONS
C A - N X (N+1) MATRIX OF SYSTEM
C B - RIGHT-HAND SIDES
C M - IF (M.EQ.0) AFTER RETURN THEN RANK(A).LT.N
C PREF(I) - PREFERENCE NUMBER FOR X(I) TO BE INDEPENDENT VARIABLE,
C 0.0.LE.PREF(I).LE.1.0, THE LOWER IS PREF(I) THE HIGHER IS
C PREFERENCE OF X(I).
C BETA(I) - COEFFICIENTS IN EXPLICIT DEPENDENCES OBTAINED IN FORM -
C          X(I)=B(I)+BETA(I)*X(K), I.NE.K.
C K - RESULTING INDEX OF INDEPENDENT VARIABLE
c
      DIMENSION A(NN,MM), B(MM), PREF(MM), BETA(MM), Y(11), X(11),
     * IRR(11), IRK(11)

      N1 = N + 1
      ID = 1
      M = 1
      DO 10 I=1,N1
        IRK(I) = 0
        IRR(I) = 0
   10 CONTINUE
   20 IR = 1
      IS = 1
      AMAX = 0.0
      DO 60 I=1,N
        IF (IRR(I)) 60, 30, 60
   30   DO 50 J=1,N1
          P = PREF(J)*ABS(A(I,J))
          IF (P-AMAX) 50, 50, 40
   40     IR = I
          IS = J
          AMAX = P
   50   CONTINUE
   60 CONTINUE
      IF (AMAX.NE.0.0) GO TO 70
c
C THIS CONDITION FOR SINGULARITY OF MATRIX MUST BE SPECIFIED
C MORE EXACTLY WITH RESPECT TO COMPUTER ACTUALLY USED
c
      M = 0
      GO TO 150
   70 IRR(IR) = IS
      DO 90 I=1,N
        IF (I.EQ.IR .OR. A(I,IS).EQ.0.0) GO TO 90
        P = A(I,IS)/A(IR,IS)
        DO 80 J=1,N1
          A(I,J) = A(I,J) - P*A(IR,J)
   80   CONTINUE
        A(I,IS) = 0.0
        B(I) = B(I) - P*B(IR)
   90 CONTINUE
      ID = ID + 1
      IF (ID.LE.N) GO TO 20
      DO 100 I=1,N
        IR = IRR(I)
        X(IR) = B(I)/A(I,IR)
        IRK(IR) = 1
  100 CONTINUE
      DO 110 K=1,N1
        IF (IRK(K).EQ.0) GO TO 120
  110 CONTINUE
  120 DO 130 I=1,N
        IR = IRR(I)
        Y(IR) = -A(I,K)/A(I,IR)
  130 CONTINUE
      DO 140 I=1,N1
        B(I) = X(I)
        BETA(I) = Y(I)
  140 CONTINUE
      B(K) = 0.0
      BETA(K) = 0.0
  150 RETURN
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
