      SUBROUTINE CDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
     8   NDE,NQ,T,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D)
C***BEGIN PROLOGUE  CDCOR
C***REFER TO  CDRIV3
C  Subroutine CDCOR is called to compute corrections to the Y array.
C  In the case of functional iteration, update Y directly from the
C  result of the last call to F.
C  In the case of the chord method, compute the corrector error and
C  solve the linear system with that as right hand side and DFDY as
C  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
C  or 5.
C***ROUTINES CALLED  CGESL,CGBSL,SCNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDCOR
      EXTERNAL FA
      COMPLEX A(MATDIM,*), DFDY(MATDIM,*), SAVE1(*), SAVE2(*), Y(*),
     8        YH(N,*), YWT(*)
      REAL D, EL(13,12), H, SCNRM2, T
      INTEGER IPVT(*)
      LOGICAL EVALFA
C***FIRST EXECUTABLE STATEMENT  CDCOR
      IF (MITER .EQ. 0) THEN
        DO 100 I = 1,N
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        D = SCNRM2(N, SAVE1, 1)/SQRT(REAL(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
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
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL CGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        DO 200 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 200      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
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
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL CGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        DO 300 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 300      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SCNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      END IF
      END
      SUBROUTINE CDCST (MAXORD,MINT,ISWFLG,EL,TQ)
C***BEGIN PROLOGUE  CDCST
C***REFER TO  CDRIV3
C  CDCST is called by CDNTL and sets coefficients used by the core
C  integrator CDSTP.  The array EL determines the basic method.
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
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDCST
      REAL EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
C***FIRST EXECUTABLE STATEMENT  CDCST
      FACTRL(1) = 1.E0
      IF (MAXORD .GE. 2) THEN
        DO 10 I = 2,MAXORD
 10       FACTRL(I) = REAL(I)*FACTRL(I-1)
      END IF
C                                             Compute Adams coefficients
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.E0
        DO 40 I = 1,MAXORD+1
          SUM = 0.E0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/REAL(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        EL(2,2) = 1.E0
        EL(3,2) = 1.E0
        IF (MAXORD .GE. 3) THEN
          DO 60 J = 3,MAXORD
            EL(2,J) = REAL(J-1)*EL(2,J-1)
            DO 50 I = 3,J
 50           EL(I,J) = REAL(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60         EL(J+1,J) = 1.E0
        END IF
        IF (MAXORD .GE. 2) THEN
          DO 80 J = 2,MAXORD
            EL(1,J) = EL(1,J-1) + GAMMA(J)
            EL(2,J) = 1.E0
            DO 80 I = 3,J+1
 80           EL(I,J) = EL(I,J)/(REAL(I-1)*FACTRL(J-1))
        END IF
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.E0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.E0/GAMMA(J+1)
 100      TQ(3,J) = -1.E0/GAMMA(J+2)
C                                              Compute Gear coefficients
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        IF (MAXORD .GE. 2) THEN
          DO 130 J = 2,MAXORD
            EL(1,J) = REAL(J)*EL(1,J-1)
            DO 120 I = 2,J
 120          EL(I,J) = REAL(J)*EL(I,J-1) + EL(I-1,J-1)
 130        EL(J+1,J) = 1.E0
          SUM = 1.E0
          DO 150 J = 2,MAXORD
            SUM = SUM + 1.E0/REAL(J)
            DO 150 I = 1,J+1
 150          EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        END IF
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.E0/FACTRL(J-1)
          TQ(2,J) = REAL(J+1)/EL(1,J)
 170      TQ(3,J) = REAL(J+2)/EL(1,J)
      END IF
C                          Compute constants used in the stiffness test.
C                          These are the ratio of TQ(2,NQ) for the Gear
C                          methods to those for the Adams methods.
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.E0
          DO 190 I = 1,MXRD
            SUM = 0.E0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/REAL(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        IF (MXRD .GE. 2) THEN
          SUM = 1.E0
          DO 200 I = 2,MXRD
            SUM = SUM + 1.E0/REAL(I)
 200        EL(1+I,1) = -REAL(I+1)*SUM*GAMMA(I+1)
        END IF
      END IF
      END
      SUBROUTINE CDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,Y,YWT,H,MNTOLD,MTROLD,NFE,RC,YH,
     8   A,CONVRG,EL,IER,IPVT,NQ,NWAIT,RH,RMAX,SAVE2,TQ,TREND,ISWFLG)
C***BEGIN PROLOGUE  CDNTL
C***REFER TO  CDRIV3
C  Subroutine CDNTL is called to set parameters on the first call
C  to CDSTP, on an internal restart, or when the user has altered
C  MINT, MITER, and/or H.
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, CDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C***ROUTINES CALLED  CDCST,CDSCL,CGEFA,CGESL,CGBFA,CGBSL,SCNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDNTL
      EXTERNAL F,FA
      COMPLEX A(MATDIM,*), SAVE1(*), SAVE2(*), Y(*), YH(N,*), YWT(*)
      REAL EL(13,12), EPS, H, HMAX, HNEW, HOLD, OLDL0, RC, RH, RMAX,
     8     RMINIT, SCNRM2, SMAX, SMIN, SUM, SUM0, T, TQ(3,12), TREND
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.E0)
C***FIRST EXECUTABLE STATEMENT  CDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) CALL CDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
        RC = 0.E0
        CONVRG = .FALSE.
        TREND = 1.E0
        RMAX = RMINIT
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              CALL CGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL CGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              CALL CGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL CGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            DO 150 I = 1,NDE
              IF(A(I,1) .EQ. CMPLX(0.E0)) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = CMPLX(0.E0)
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = SCNRM2(NDE, SAVE1, 1)
        SUM0 = 1.E0/MAX(1.E0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.E0 + (SMIN/SMAX)**2)/SQRT(REAL(NDE))
        H = SIGN(MIN(2.E0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.E0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL CDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          HNEW = H
          RH = H/HOLD
          H = HOLD
          CALL CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          H = SIGN(MIN(ABS(H), ABS(HNEW)), H)
        END IF
      END IF
      END
      SUBROUTINE CDNTP (H,K,N,NQ,T,TOUT,YH,Y)
C***BEGIN PROLOGUE  CDNTP
C***REFER TO  CDRIV3
C   Subroutine CDNTP interpolates the K-th derivative of Y at TOUT,
C   using the data in the YH array.  If K has a value greater than NQ,
C   the NQ-th derivative is calculated.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDNTP
      COMPLEX Y(*), YH(N,*)
      REAL FACTOR, H, R, T, TOUT
C***FIRST EXECUTABLE STATEMENT  CDNTP
      KUSED = MIN(K, NQ)
      IF (KUSED .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        FACTOR = 1.E0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*REAL(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        IF (KUSED .NE. NQ) THEN
          R = ((TOUT - T)/H)
          DO 80 JJ = KUSED+1,NQ
            J = K + 1 + NQ - JJ
            FACTOR = 1.E0
            DO 60 KK = 1,KUSED
 60           FACTOR = FACTOR*REAL(J-KK)
            DO 70 I = 1,N
 70           Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80         CONTINUE
        END IF
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
      SUBROUTINE CDPSC (KSGN,N,NQ,YH)
C***BEGIN PROLOGUE  CDPSC
C***REFER TO  CDRIV3
C     This subroutine computes the predicted YH values by effectively
C     multiplying the YH array by the Pascal triangle matrix when KSGN
C     is +1, and performs the inverse function when KSGN is -1.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDPSC
      COMPLEX YH(N,*)
C***FIRST EXECUTABLE STATEMENT  CDPSC
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
      SUBROUTINE CDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
     8   NQ,SAVE2,T,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,IER,IPVT,SAVE1,
     8   ISWFLG,BND)
C***BEGIN PROLOGUE  CDPST
C***REFER TO  CDRIV3
C  Subroutine CDPST is called to reevaluate the partials.
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C***ROUTINES CALLED  CGEFA,CGBFA,SCNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDPST
      EXTERNAL F,FA,JACOBN
      COMPLEX A(MATDIM,*), CFCTR, DFDY(MATDIM,*), DY, SAVE1(*),
     8        SAVE2(*), Y(*), YH(N,*), YJ, YWT(*)
      REAL BND, DFDYMX, EL(13,12), EPSJ, ETA, ETATST, H, RFCTR,
     8     SCNRM2, T, UROUND, ZMAX, ZMIN
      INTEGER IPVT(*)
      LOGICAL IER, LOOP
       PARAMETER(ETATST = .5E0, ITERMX = 3)
C***FIRST EXECUTABLE STATEMENT  CDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (ISWFLG .EQ. 3) BND = SCNRM2(N*N, DFDY, 1)
          RFCTR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = RFCTR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          EPSJ = UROUND**(1.E0/3.E0)
          DO 170 J = 1,N
            IF (ABS(Y(J)).GT.ABS(YWT(J))) THEN
              DY = EPSJ*Y(J)
            ELSE
              DY = EPSJ*YWT(J)
            END IF
            IF (DY .EQ. CMPLX(0.E0)) DY = CMPLX(EPSJ, EPSJ)
            ITER = 0
 120        YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            Y(J) = YJ
            NFE = NFE + 1
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              DO 130 I = 1,N
                IF (SAVE1(I) .NE. SAVE2(I)) THEN
                  ETA = ABS(SAVE2(I))*UROUND/
     8            (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                  IF (ETA .GE. ETATST) THEN
                    DY = DY*10.E0
                    GO TO 120
                  END IF
                END IF
 130            CONTINUE
            END IF
            CFCTR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*CFCTR
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = SCNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL CGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          RFCTR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = RFCTR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          EPSJ = UROUND**(1.E0/3.E0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 265 K = J,N,MW
              IF (ABS(Y(K)).GT.ABS(YWT(K))) THEN
                DY = EPSJ*Y(K)
              ELSE
                DY = EPSJ*YWT(K)
              END IF
              IF (DY .EQ. CMPLX(0.E0)) DY = CMPLX(EPSJ, EPSJ)
              DFDY(MW,K) = Y(K)
 265          Y(K) = Y(K) + DY
            ITER = 0
 270        CALL F (N, T, Y, SAVE1)
            NFE = NFE + 0
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              LOOP = .FALSE.
              DO 290 K = J,N,MW
                I1 = MAX(1, K-MU)
                I2 = MIN(K+ML, N)
                DO 280 I = I1,I2
                  IF (SAVE1(I) .NE. SAVE2(I)) THEN
                    ETA = ABS(SAVE2(I))*UROUND/
     8              (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                    IF (ETA .GE. ETATST) THEN
                      DY = (Y(K) - DFDY(MW,K))*10.E0
                      Y(K) = DFDY(MW,K) + DY
                      LOOP = .TRUE.
                      GO TO 290
                    END IF
                  END IF
 280              CONTINUE
 290            CONTINUE
              IF (LOOP) GO TO 270
            END IF
            DO 330 K = J,N,MW
              DY = Y(K) - DFDY(MW,K)
              Y(K) = DFDY(MW,K)
              CFCTR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = CFCTR*(SAVE1(I3) - SAVE2(I3))
 330          CONTINUE
 340        CONTINUE
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.E0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 345 I = I1,I2
              ZMAX = MAX(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
              ZMIN = MIN(ABS(REAL(DFDY(I,J))), ABS(AIMAG(DFDY(I,J))))
              IF (ZMAX .NE. 0.E0)
     8        DFDYMX = MAX(DFDYMX, ZMAX*SQRT(1.E0+ (ZMIN/ZMAX)**2))
 345          CONTINUE
          BND = 0.E0
          IF (DFDYMX .NE. 0.E0) THEN
            DO 350 J = 1,N
              I1 = MAX(ML+1, MW+1-J)
              I2 = MIN(MW+N-J, MW+ML)
              DO 350 I = I1,I2
                BND = BND + (REAL(DFDY(I,J))/DFDYMX)**2 +
     8          (AIMAG(DFDY(I,J))/DFDYMX)**2
 350            CONTINUE
            BND = DFDYMX*SQRT(BND)/(-EL(1,NQ)*H)
          END IF
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 360 J = 1,N
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL CGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
      END IF
      END
      SUBROUTINE CDRIV1 (N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)
C***BEGIN PROLOGUE  CDRIV1
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             COMPLEX VALUED
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of CDRIV1 is to solve N (200 or fewer)
C            ordinary differential equations of the form
C            dY(I)/dT = F(Y(I),T), given the initial conditions
C            Y(I) = YI.  CDRIV1 allows complex-valued differential
C            equations.
C***DESCRIPTION
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
C           real, independent variable.
C
C    As an aid in selecting the proper program, the following is a
C    discussion of the important options or restrictions associated with
C    each program:
C
C      A. CDRIV1 should be tried first for those routine problems with
C         no more than 200 differential equations.  Internally this
C         routine has two important technical defaults:
C           1. Numerical approximation of the Jacobian matrix of the
C              right hand side is used.
C           2. The stiff solver option is used.
C         Most users of CDRIV1 should not have to concern themselves
C         with these details.
C
C      B. CDRIV2 should be considered for those problems for which
C         CDRIV1 is inadequate (CDRIV2 has no explicit restriction on
C         the number of differential equations.)  For example, CDRIV1
C         may have difficulty with problems having zero initial
C         conditions and zero derivatives.  In this case CDRIV2, with an
C         appropriate value of the parameter EWT, should perform more
C         efficiently.  CDRIV2 provides three important additional
C         options:
C           1. The nonstiff equation solver (as well as the stiff
C              solver) is available.
C           2. The root-finding option is available.
C           3. The program can dynamically select either the non-stiff
C              or the stiff methods.
C         Internally this routine also defaults to the numerical
C         approximation of the Jacobian matrix of the right hand side.
C
C      C. CDRIV3 is the most flexible, and hence the most complex, of
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
C    The function of CDRIV1 is to solve N (200 or fewer) ordinary
C    differential equations of the form dY(I)/dT = F(Y(I),T), given the
C    initial conditions Y(I) = YI.  CDRIV1 is to be called once for each
C    output point.
C
C  III.  PARAMETERS  ...................................................
C
C    The user should use parameter names in the call sequence of CDRIV1
C    for those quantities whose value may be altered by CDRIV1.  The
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
C             in the call sequence of the user-provided routine F.
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
C                  user.  Unless CDRIV1 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, CDRIV1 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling CDRIV1
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling CDRIV1
C                  again.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  On output, the adjusted relative accuracy if
C             the input value was too small.  The value of EPS should be
C             set as large as is reasonable, because the amount of work
C             done by CDRIV1 increases as EPS decreases.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW complex words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       COMPLEX WORK(...)
C             The length of WORK should be at least N*N + 10*N + 225
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to CDRIV1.
C
C***LONG DESCRIPTION
C
C  IV.  USAGE  .........................................................
C
C                   PROGRAM SAMPLE
C                   COMPLEX Y(...), WORK(...)
C                   OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C                   N = ...                       Number of equations
C                   T = ...                       Initial point
C                   DO 10 I = 1,N
C              10     Y(I) = ...                  Set initial conditions
C                   TOUT = T
C                   MSTATE = 1
C                   EPS = ...
C                   LENW = ...
C              20   CALL CDRIV1 (N, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C                   IF (MSTATE .GT. 2) STOP
C                   WRITE(6, 100) TOUT, (Y(I), I=1,N)
C                   TOUT = TOUT + 1.
C                   IF (TOUT .LE. 10.) GO TO 20
C              100  FORMAT(...)
C                   END (Sample)
C
C    The user must write a subroutine called F to evaluate the right
C    hand side of the differential equations.  It is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   COMPLEX Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C    This computes YDOT = F(Y,T), the right hand side of the
C    differential equations.  Here Y is a vector of length at least N.
C    The actual length of Y is determined by the user's declaration in
C    the program which calls CDRIV1.  Thus the dimensioning of Y in F,
C    while required by FORTRAN convention, does not actually allocate
C    any storage.  When this subroutine is called, the first N
C    components of Y are intermediate approximations to the solution
C    components.  The user should not alter these values.  Here YDOT is
C    a vector of length N.  The user should only compute YDOT(I) for I
C    from 1 to N.
C
C  V.  OTHER COMMUNICATION TO THE USER  ................................
C
C    The solver communicates to the user through the parameters above.
C    In addition it writes diagnostic messages through the standard
C    error handling program XERROR.  That program will terminate the
C    user's run if it detects a probable problem setup error, e.g.,
C    insufficient storage allocated by the user for the WORK array.  For
C    further information see section III-A of the writeup for CDRIV3.
C
C  VI.  REMARKS  .......................................................
C
C    A. There are user-written routines which are only required by
C       CDRIV2 or CDRIV3 when certain parameters are set.  Thus a
C       message warning of unsatisfied externals may be issued during
C       the load or link phase.  This message can be ignored unless it
C       refers to F.
C
C    For other information, see section IV of the writeup for CDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  CDRIV3,R1MACH,XERROR
C***END PROLOGUE  CDRIV1
      EXTERNAL F, JACOBN, FA, G
      COMPLEX WORK(*), Y(*)
      REAL EPS, EWT, G, HMAX, R1MACH, T, TOUT
      PARAMETER(MXN = 200, IDLIW = 21, MXLIW = IDLIW + MXN)
      INTEGER IWORK(MXLIW)
      CHARACTER MSG*103
      PARAMETER(NROOT = 0, EWT = 1.E0, IERROR = 2, MINT = 2, MITER = 2,
     8          IMPL = 0, ML = 0, MU = 0, MXORD = 5, NDE = 0,
     8          MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  CDRIV1
      IF (N .GT. MXN) THEN
        WRITE(MSG, '(''CDRIV115FE Illegal input.  The number of '',
     8  ''equations,'', I8, '', is greater than the maximum allowed.'')
     8  ') N
        CALL XERROR(MSG, 97, 15, 2)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = SQRT(R1MACH(2))
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 9*N + 204)) THEN
        LNWCHK = N*N + 9*N + 204 + LENIW
        WRITE(MSG, '(''CDRIV116FE Insufficient storage allocated for '',
     8  ''the work array.  The required storage is at least'', I8)')
     8  LNWCHK
        CALL XERROR(MSG, 103, 16, 2)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
          II = I + LENWCM
 20       IWORK(I) = INT(REAL(WORK(II)))
      END IF
      CALL CDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWT,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      DO 40 I = 1,LENIW
        II = LENWCM + I
 40     WORK(II) = CMPLX(REAL(IWORK(I)))
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE CDRIV2 (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     8   LENW,IWORK,LENIW,G)
C***BEGIN PROLOGUE  CDRIV2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             COMPLEX VALUED
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of CDRIV2 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  CDRIV2 allows complex-valued differential
C            equations.
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The function of CDRIV2 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.
C    CDRIV2 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of CDRIV2
C    for those quantities whose value may be altered by CDRIV2.  The
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
C             G.
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   COMPLEX Y(*), YDOT(*)
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
C             user's declaration in the program which calls CDRIV2.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.
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
C                  user.  Unless CDRIV2 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, CDRIV2 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling CDRIV2
C                  again.  Other than an error in problem setup, the
C                  most likely cause for this condition is trying to
C                  integrate a stiff set of equations with the non-stiff
C                  integrator option. (See description of MINT below.)
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling CDRIV2
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling CDRIV2 again.
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  CDRIV2 attempts
C             to find the value of T at which one of the equations
C             changes sign.  CDRIV2 can find at most one root per
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
C             because the amount of work done by CDRIV2 increases as
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
C             WORK is an array of LENW complex words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       COMPLEX WORK(...)
C             The length of WORK should be at least
C               16*N + 2*NROOT + 204         if MINT is 1, or
C               N*N + 9*N + 2*NROOT + 204    if MINT is 2, or
C               N*N + 16*N + 2*NROOT + 204   if MINT is 3,
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to CDRIV2.
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
C             of IWORK should not be disturbed between calls to CDRIV2.
C
C    G      = A real FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to
C             obtain the value of each of the NROOT equations for which
C             a root is desired.  G is of the form:
C                   REAL FUNCTION G (N, T, Y, IROOT)
C                   COMPLEX Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls CDRIV2.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.
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
C         NSTEP    The number of steps taken.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. On any return from CDRIV2 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    B. There are user-written routines which are only required by
C       CDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if it refers to G and NROOT is 0.  A reference to any
C       other unsatisfied external can be ignored.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to CDRIV2.
C
C  V.  USAGE  ..........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               COMPLEX WORK(...), Y(...)        See II. for
C               INTEGER IWORK(...)               required dimensions for
C                                                WORK and IWORK
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C               N = ...                          Number of equations
C               T = ...                          Initial point
C               DO 10 I = 1,N
C          10     Y(I) = ...                     Set initial conditions
C               TOUT = T
C               MSTATE = 1
C               NROOT = 0
C               EPS = ...
C               EWT = ...
C               MINT = 1
C               LENW = ...
C               LENIW = ...
C          20   CALL CDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, G)
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  CDRIV3,R1MACH,XERROR
C***END PROLOGUE  CDRIV2
      EXTERNAL F, JACOBN, FA, G
      COMPLEX WORK(*), Y(*)
      REAL EPS, EWT, EWTCOM(1), G, HMAX, R1MACH, T, TOUT
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, ML = 0, MU = 0, NDE = 0, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  CDRIV2
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(MSG, '(''CDRIV21FE Illegal input.  Improper value for '',
     8  ''the integration method flag,'', I8)') MINT
        CALL XERROR(MSG, 81, 21, 2)
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
      IF (EWT .NE. 0.E0) THEN
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
      HMAX = SQRT(R1MACH(2))
      CALL CDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE CDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G)
C***BEGIN PROLOGUE  CDRIV3
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             COMPLEX VALUED
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of CDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  CDRIV3
C            allows complex-valued differential equations.
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of CDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, CDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    CDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of CDRIV3
C    for those quantities whose value may be altered by CDRIV3.  The
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
C             JACOBN, FA, USERS, and G.
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   COMPLEX Y(*), YDOT(*)
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
C             user's declaration in the program which calls CDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, CDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling CDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling CDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling CDRIV3 again.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means CDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means CDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means CDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  CDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  CDRIV3 can find at most one root per
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
C             because the amount of work done by CDRIV3 increases as EPS
C             decreases.
C
C    EWT    = (Input) Problem zero, i.e., the smallest, nonzero,
C             physically meaningful value for the solution.  (Array,
C             possibly of length one.  See following description of
C             IERROR.)
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
C                          by the user-written routine named USERS.
C                          This option allows all matrix algebra and
C                          storage decisions to be made by the user.
C                          The routine USERS is called by CDRIV3 when
C                          certain linear systems must be solved.  The
C                          user may choose any method to form, store and
C                          solve these systems in order to obtain the
C                          solution result that is returned to CDRIV3.
C                          In particular, this allows sparse matrix
C                          methods to be used.
C                          The call sequence for this routine is
C
C                           SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2,
C                          8              T, H, EL, IMPL, N, NDE, IFLAG)
C                           COMPLEX Y(*), YH(*), YWT(*),
C                          8        SAVE1(*), SAVE2(*)
C
C                          The input variable IFLAG indicates what
C                          action is to be taken.  Subroutine USERS
C                          should perform the following operations,
C                          depending on the value of IFLAG and IMPL.
C
C                          IFLAG = 0
C                            IMPL = 0.  USERS is not called.
C                            IMPL = 1 or 2.  Solve the system
C                                A*X = SAVE2,
C                              returning the result in SAVE2.  The array
C                              SAVE1 can be used as a work array.
C
C                          IFLAG = 1
C                            IMPL = 0.  Compute, decompose and store the
C                              matrix (I - H*EL*J), where I is the
C                              identity matrix and J is the Jacobian
C                              matrix of the right hand side.  The array
C                              SAVE1 can be used as a work array.
C                            IMPL = 1 or 2. Compute, decompose and store
C                              the matrix (A - H*EL*J).  The array SAVE1
C                              can be used as a work array.
C
C                          IFLAG = 2
C                            IMPL = 0.   Solve the system
C                                (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                              returning the result in SAVE2.
C                            IMPL = 1, or 2.  Solve the system
C                              (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                              returning the result in SAVE2.
C                            The array SAVE1 should not be altered.
C
C                          When using a value of MITER = 3, the
C                          subroutine FA is not required, even if IMPL
C                          is not 0.  For further information on using
C                          this option, see section IV-F below.
C
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
C             WORK is an array of LENW complex words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       COMPLEX WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             CDRIV3.
C
C      IMPL =   0                   1                   2
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N +       Not allowed         Not allowed
C              2*NROOT + 204
C
C         1,2  N*N+(MXORD+4)*N     2*N*N+(MXORD+4)*N   N*N+(MXORD+5)*N
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204
C
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +
C              (MXORD+5)*N +       (MXORD+6)*N +       (MXORD+6)*N +
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
C             of IWORK should not be disturbed between calls to CDRIV3.
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
C             approximated numerically by CDRIV3, and we only ask the
C             user to tell CDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   COMPLEX Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls CDRIV3.  Thus the dimensioning of Y in
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
C             K=I-J+MU+1.
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
C                   COMPLEX Y(*), A(MATDIM,*)
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
C                   COMPLEX Y(*), A(*)
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
C             program which calls CDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to CDRIV3.
C
C    G      = A real FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   REAL FUNCTION G (N, T, Y, IROOT)
C                   COMPLEX Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls CDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.
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
C       come from CDRIV3:
C
C        No.  Type         Message
C        ---  ----         -------
C         1   Fatal        From CDRIV2: The integration method flag has
C                          an illegal value.
C         2   Warning      The output point is inconsistent with the
C                          value of NTASK and T.
C         3   Warning      Number of steps to reach TOUT exceeds MXSTEP.
C         4   Recoverable  Requested accuracy is too stringent.
C         5   Warning      Step size is below the roundoff level.
C         6   Fatal        EPS is less than zero.
C         7   Fatal        N is not positive.
C         8   Fatal        Insufficient work space provided.
C         9   Fatal        Improper value for MINT, MITER and/or IMPL.
C        10   Fatal        The IWORK array is too small.
C        11   Fatal        The step size has gone to zero.
C        12   Fatal        Excessive amount of work.
C        13   Fatal        For IMPL=1 or 2, the matrix A is singular.
C        14   Fatal        MXORD is not positive.
C        15   Fatal        From CDRIV1: N is greater than 200.
C        16   Fatal        From CDRIV1: The WORK array is too small.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         CDNTP, CDZRO, CDSTP, CDNTL, CDPST, CDCOR, CDCST,
C         CDPSC, and CDSCL;
C         CGEFA, CGESL, CGBFA, CGBSL, and SCNRM2 (from LINPACK)
C         R1MACH (from the Bell Laboratories Machine Constants Package)
C         XERROR (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from CDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. There are user-written routines which are only required by
C       CDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if: it refers to JACOBN and MITER is not 1 or 4, or
C       it refers to FA and IMPL is 0 or MITER is 3, or it refers to
C       USERS and MITER is not 3, or it refers to G and NROOT is 0.
C
C    D. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to CDRIV3.
C
C    E. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to CDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to CDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    F. As the price for complete control of matrix algebra, the CDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C             REAL EPSJ, H, R1MACH, T, UROUND
C             COMPLEX DFDY(N,N), R, SAVE1(N), SAVE2(N), Y(N), YJ, YWT(N)
C             UROUND = R1MACH(4)
C             EPSJ = UROUND**(1.E0/3.E0)
C             DO 30 J = J1,J2
C               IF (ABS(Y(J)).GT.ABS(YWT(J))) THEN
C                 R = EPSJ*Y(J)
C               ELSE
C                 R = EPSJ*YWT(J)
C               END IF
C               IF (R .EQ. CMPLX(0.E0)) R = CMPLX(EPSJ, EPSJ)
C               YJ = Y(J)
C               Y(J) = Y(J) + R
C               CALL F (N, T, Y, SAVE1)
C               Y(J) = YJ
C               DO 20 I = I1,I2
C        20       DFDY(I,J) = (SAVE1(I) - SAVE2(I))/R
C        30     CONTINUE
C       Many problems give rise to structured sparse Jacobians, e.g.,
C       block banded.  It is possible to approximate them with fewer
C       function evaluations than the above procedure uses; see Curtis,
C       Powell and Reid, J. Inst. Maths Applics, (1974), Vol. 13,
C       pp. 117-119.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  CDSTP,CDNTP,CDZRO,CGEFA,CGESL,CGBFA,CGBSL,
C                    SCNRM2,R1MACH,XERROR
C***END PROLOGUE  CDRIV3
      EXTERNAL F, JACOBN, FA, G
      COMPLEX WORK(*), Y(*)
      REAL AE, AVGH, AVGORD, BIG, EL(13,12), EPS, EWT(*), G,
     8     GLAST, GNOW, H, HMAX, HOLD, HSIGN, HUSED, NROUND, RC, RE,
     8     RMAX, R1MACH, SCNRM2, SIZE, SUM, T, TLAST, TOUT, TQ(3,12),
     8     TREND, TROOT, UROUND
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.E0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,
     8          INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20)
      PARAMETER(INDPRT = 21, INDPVT = 22)
C***FIRST EXECUTABLE STATEMENT  CDRIV3
      UROUND = R1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = R1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.E0) THEN
        WRITE(MSG, '(''CDRIV36FE Illegal input.  EPS,'', E16.8,
     8  '', is negative.'')') EPS
        CALL XERROR(MSG, 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''CDRIV37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
        CALL XERROR(MSG, 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''CDRIV314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
        CALL XERROR(MSG, 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1)) THEN
        WRITE(MSG, '(''CDRIV39FE Illegal input.  Improper value for '',
     8  ''MINT, MITER and/or IMPL.'')')
        CALL XERROR(MSG, 69, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(MSG, '(''CDRIV310FE Illegal input.  Insufficient '',
     8  ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LIWCHK
        CALL XERROR(MSG, 164, 10, 2)
        RETURN
      END IF
C                                                Allocate the work array
C                                         IYH is the index of YH in WORK
      IF (MINT .EQ. 1 .OR. MINT .EQ. 3) THEN
        MAXORD = MIN(MXORD, 12)
      ELSE IF (MINT .EQ. 2) THEN
        MAXORD = MIN(MXORD, 5)
      END IF
      IDFDY = IYH + (MAXORD + 1)*N
C                                             IDFDY is the index of DFDY
C
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        IYWT = IDFDY
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IYWT = IDFDY + N*N
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IYWT = IDFDY + (2*ML + MU + 1)*N
      END IF
C
C                                               IYWT is the index of YWT
      ISAVE1 = IYWT + N
C                                           ISAVE1 is the index of SAVE1
      ISAVE2 = ISAVE1 + N
C                                           ISAVE2 is the index of SAVE2
      IGNOW = ISAVE2 + N
C                                             IGNOW is the index of GNOW
      ITROOT = IGNOW + NROOT
C                                           ITROOT is the index of TROOT
      IA = ITROOT + NROOT
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
        WRITE(MSG, '(''CDRIV38FE Illegal input.  Insufficient '',
     8  ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LENCHK
        CALL XERROR(MSG, 162, 8, 2)
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
        IF (T .EQ. TOUT) RETURN
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
        WORK(IHMAX) = CMPLX(HMAX)
        H = (TOUT - T)*(1.E0 - 4.E0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = CMPLX(H)
        HSIGN = SIGN(1.E0, H)
        WORK(IHSIGN) = CMPLX(HSIGN)
        IWORK(IJTASK) = 0
        AVGH = 0.E0
        AVGORD = 0.E0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        WORK(IT) = CMPLX(T)
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
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
      AVGH = REAL(WORK(IAVGH))
      AVGORD = REAL(WORK(IAVGRD))
      HOLD = REAL(WORK(IHOLD))
      RC = REAL(WORK(IRC))
      RMAX = REAL(WORK(IRMAX))
      TREND = REAL(WORK(ITREND))

      DO J = 1,12
        JEL1 = IEL + (J-1)*13 - 1
        DO I = 1,13
          JEL = JEL1 + I
          EL(I,J) = REAL(WORK(JEL))
        end do
      end do

      DO J = 1,12
        JTQ1 = ITQ + (J-1)*3 - 1
        DO I = 1,3
          JTQ = JTQ1 + I
          TQ(I,J) = REAL(WORK(JTQ))
        end do
      end do
      T = REAL(WORK(IT))
      H = REAL(WORK(IH))
      HSIGN = REAL(WORK(IHSIGN))
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
      IF(NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. REAL(WORK(ITOUT))*HSIGN) THEN
              TROOT = REAL(WORK(ITOUT))
              CALL CDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
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
              IF (REAL(WORK(JTROOT))*HSIGN .LE. TROOT*HSIGN) THEN
C
C                                              Check for multiple roots.
c
                IF (REAL(WORK(JTROOT)) .EQ. REAL(WORK(ITOUT)) .AND.
     8          I .GT. IWORK(INROOT)) THEN
                  IROOT = I
                  TROOT = REAL(WORK(JTROOT))
                  GO TO 60
                END IF
                IF (REAL(WORK(JTROOT))*HSIGN .GT.
     8          REAL(WORK(ITOUT))*HSIGN) THEN
                  IROOT = I
                  TROOT = REAL(WORK(JTROOT))
                END IF
              END IF
 50           CONTINUE
 60         IWORK(INROOT) = IROOT
            WORK(ITOUT) = CMPLX(TROOT)
            IWORK(IJROOT) = NTASK
            IF (NTASK .EQ. 1) THEN
              IF (IROOT .EQ. 0) THEN
                IWORK(IJROOT) = 0
              ELSE
                IF (TOUT*HSIGN .GE. TROOT*HSIGN) THEN
                  CALL CDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
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
                CALL CDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
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
          CALL CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''CDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
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
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = CMPLX(H)
          IF (H .EQ. 0.E0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''CDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = CMPLX(H)
          IF (H .EQ. 0.E0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      END IF
C                         Implement changes in MINT, MITER, and/or HMAX.
C
      IF ((MINT .NE. IWORK(IMNTLD) .OR. MITER .NE. IWORK(IMTRLD)) .AND.
     8  MINT .NE. 3 .AND. IWORK(IMNTLD) .NE. 3) IWORK(IJTASK) = -1
      IF (HMAX .NE. REAL(WORK(IHMAX))) THEN
        H = SIGN(MIN(ABS(H), HMAX), H)
        IF (H .NE. REAL(WORK(IH))) THEN
          IWORK(IJTASK) = -1
          WORK(IH) = CMPLX(H)
        END IF
        WORK(IHMAX) = CMPLX(HMAX)
      END IF
 180  NSTEPL = IWORK(INSTEP)
      DO 190 I = 1,N
        JYH = IYH + I - 1
 190    Y(I) = WORK(JYH)
      IF (NROOT .NE. 0) THEN
        DO 200 I = 1,NROOT
          JGNOW = IGNOW + I - 1
 200      WORK(JGNOW) = CMPLX(G (N, T, Y, I))
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = CMPLX(1.E0)
        GO TO 410
      ELSE IF (IERROR .EQ. 5) THEN
        DO 250 I = 1,N
          JYWT = I + IYWT - 1
 250      WORK(JYWT) = CMPLX(EWT(I))
        GO TO 410
      END IF
C                                       Reset YWT array.  Looping point.
 260  IF (IERROR .EQ. 2) THEN
        DO 280 I = 1,N
          IF (Y(I) .EQ. CMPLX(0.E0)) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = Y(I)
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (N, T, Y, WORK(ISAVE2))
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y, WORK(IYH), WORK(IYWT), WORK(ISAVE1),
     8                 WORK(ISAVE2), T, H, WORK(IEL), IMPL, N, NDECOM,
     8                 IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              CALL CGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL CGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (N, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              CALL CGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL CGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF(WORK(JA) .EQ. CMPLX(0.E0)) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. CMPLX(0.E0)) THEN
            WORK(JYWT) = Y(J)
          ELSE
            IF (IWORK(IJTASK) .EQ. 0) THEN
              JSAVE2 = J + ISAVE2 - 1
              WORK(JYWT) = H*WORK(JSAVE2)
            ELSE
              JHYP = J + IYH + N - 1
              WORK(JYWT) = WORK(JHYP)
            END IF
          END IF
          IF (WORK(JYWT) .EQ. CMPLX(0.E0))
     8    WORK(JYWT) = CMPLX(UROUND)
 360      CONTINUE
      ELSE IF (IERROR .EQ. 3) THEN
        DO 380 I = 1,N
          JYWT = I + IYWT - 1
 380      WORK(JYWT) = CMPLX(MAX(EWT(1), ABS(Y(I))))
      ELSE IF (IERROR .EQ. 4) THEN
        DO 400 I = 1,N
          JYWT = I + IYWT - 1
 400      WORK(JYWT) = CMPLX(MAX(EWT(I), ABS(Y(I))))
      END IF
C
 410  DO 420 I = 1,N
        JYWT = I + IYWT - 1
        JSAVE2 = I + ISAVE2 - 1
 420    WORK(JSAVE2) = Y(I)/WORK(JYWT)
      SUM = SCNRM2(N, WORK(ISAVE2), 1)/SQRT(REAL(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.E0 + 10.E0*UROUND)
        WRITE(MSG, '(''CDRIV34REC At T,'', E16.8, '', the requested '',
     8  ''accuracy, EPS, was not obtainable with the machine '',
     8  ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(E16.8)') EPS
        CALL XERROR(MSG, 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(MSG, '(''CDRIV35WRN At T,'', E16.8, '', the step size,'',
     8  E16.8, '', is smaller than the roundoff level of T.  '')') T, H
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
     8  ''change in the right hand side of the differential '',
     8  ''equations.'')')
        CALL XERROR(MSG, 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''CDRIV33WRN At T,'', E16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', E16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG, 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL CDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL CDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, N,
     8            NDECOM, WORK(IYWT), UROUND, AVGH, AVGORD, H, HUSED,
     8            IWORK(IJTASK), IWORK(IMNTLD), IWORK(IMTRLD),
     8            IWORK(INFE), IWORK(INJE), IWORK(INQUSD),
     8            IWORK(INSTEP), T, Y, WORK(IYH), WORK(IA), CONVRG,
     8            WORK(IDFDY), EL, HOLD, IWORK(INDPVT), JSTATE,
     8            IWORK(INQ), IWORK(INWAIT), RC, RMAX, WORK(ISAVE1),
     8            WORK(ISAVE2), TQ, TREND, MINT, IWORK(IMTRSV),
     8            IWORK(IMXRDS))
C
      WORK(IH) = CMPLX(H)
      WORK(IT) = CMPLX(T)
      GO TO (470, 670, 680, 690), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (ABS(H) .GE. UROUND*ABS(T) .AND. NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = REAL(WORK(JGNOW))
          GNOW = G (N, T, Y, I)
          WORK(JGNOW) = CMPLX(GNOW)
          IF (GLAST*GNOW .GT. 0.E0) THEN
            WORK(JTROOT) = CMPLX(T + H)
          ELSE
            IF (GNOW .EQ. 0.E0) THEN
              WORK(JTROOT) = CMPLX(T)
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.E0) THEN
                WORK(JTROOT) = CMPLX(T + H)
              ELSE
                TLAST = T - HUSED
                IROOT = I
                TROOT = T
                CALL CDZRO (AE, G, H, N, IWORK(INQ), IROOT, RE, T,
     8                       WORK(IYH), UROUND,  TROOT, TLAST, GNOW,
     8                      GLAST,  WORK(ISAVE1))
                WORK(JTROOT) = CMPLX(TROOT)
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
            IF (REAL(WORK(JTROOT))*HSIGN .LT. TROOT*HSIGN) THEN
              TROOT = REAL(WORK(JTROOT))
              IROOT = I
            END IF
 510        CONTINUE
          IWORK(INROOT) = IROOT
          WORK(ITOUT) = CMPLX(TROOT)
          IF (TROOT*HSIGN .LE. TOUT*HSIGN) THEN
            CALL CDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
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
        CALL CDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
        T = TOUT
        GO TO 580
C                               TOUT is assumed to have been attained
C                               exactly if T is within twenty roundoff
C                               units of TOUT, relative to MAX(TOUT, T).
C
      ELSE IF (NTASK .EQ. 2) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = CMPLX(H)
            IF (H .EQ. 0.E0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = CMPLX(H)
            IF (H .EQ. 0.E0) GO TO 670
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
      WORK(IAVGH) = CMPLX(AVGH)
      WORK(IAVGRD) = CMPLX(AVGORD)
      WORK(IHUSED) = CMPLX(HUSED)
      WORK(IHOLD) = CMPLX(HOLD)
      WORK(IRC) = CMPLX(RC)
      WORK(IRMAX) = CMPLX(RMAX)
      WORK(ITREND) = CMPLX(TREND)
      DO 582 J = 1,12
        JEL1 = IEL + (J-1)*13 - 1
        DO 582 I = 1,13
          JEL = JEL1 + I
 582      WORK(JEL) = CMPLX(EL(I,J))
      DO 584 J = 1,12
        JTQ1 = ITQ + (J-1)*3 - 1
        DO 584 I = 1,3
          JTQ = JTQ1 + I
 584      WORK(JTQ) = CMPLX(TQ(I,J))
      IF (IWORK(IJTASK) .EQ. 0) RETURN
      BIG = 0.E0
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
C                                        Fatal errors are processed here
C
 670  WRITE(MSG, '(''CDRIV311FE At T,'', E16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL XERROR(MSG, 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''CDRIV312FE At T,'', E16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG, 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''CDRIV313FE At T,'', E16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG, 74, 13, 2)
      RETURN
      END
      SUBROUTINE CDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)
C***BEGIN PROLOGUE  CDSCL
C***REFER TO  CDRIV3
C   This subroutine rescales the YH array whenever the step size
C   is changed.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850319   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDSCL
      COMPLEX YH(N,*)
      REAL H, HMAX, RC, RH, RMAX, R1
C***FIRST EXECUTABLE STATEMENT  CDSCL
      IF (H .LT. 1.E0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.E0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      END
      SUBROUTINE CDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,
     8   MITER,ML,MU,N,NDE,YWT,UROUND,AVGH,AVGORD,H,HUSED,JTASK,MNTOLD,
     8   MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,HOLD,IPVT,
     8   JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,MTRSV,
     8   MXRDSV)
C***BEGIN PROLOGUE  CDSTP
C***REFER TO  CDRIV3
C
C  CDSTP performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Communication with CDSTP is done with the following variables:
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
C              the value used in CDSTP.
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
C
C***ROUTINES CALLED  CDNTL,CDPST,CDCOR,CDPSC,CDSCL,SCNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860513   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDSTP
      EXTERNAL F, JACOBN, FA
      COMPLEX A(MATDIM,*), DFDY(MATDIM,*), SAVE1(*), SAVE2(*), Y(*),
     8        YH(N,*), YWT(*)
      REAL AVGH, AVGORD, BIAS1, BIAS2, BIAS3, BND, CTEST, D, DENOM, D1,
     8     EL(13,12), EPS, ERDN, ERUP, ETEST, H, HMAX, HN, HOLD, HS,
     8     HUSED, NUMER, RC, RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL,
     8     RMNORM, SCNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, UROUND,
     8     Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3E0, BIAS2 = 1.2E0, BIAS3 = 1.4E0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3E0, RMFAIL = 2.E0,
     8          RMNORM = 10.E0, TRSHLD = 1.E0)
C***FIRST EXECUTABLE STATEMENT  CDSTP
      BND = 0.E0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL CDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
      END IF
      IER = .FALSE.
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL CDPSC (1, N, NQ,  YH)
      EVALJC = ((ABS(RC - 1.E0) .GT. RCTEST) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL CDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, IER, IPVT, SAVE1, ISWFLG, BND)
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.E0
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = CMPLX(0.E0)
C
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
 130  CALL CDCOR (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,
     8            MU, N, NDE, NQ, T, Y, YH, YWT,  EVALFA, SAVE1,
     8            SAVE2,  A, D)
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = SCNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = SCNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = SCNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.E0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.E0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9E0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.E0*TREND, 1.E0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
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
      CALL CDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3E0
      ELSE
        RH = .9E0*(EPS/CTEST)**(.2E0)
      END IF
      IF (RH*H .EQ. 0.E0) GO TO 400
      CALL CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C                          The corrector has converged.  CONVRG is set
C                          to .TRUE. if partial derivatives were used,
C                          to indicate that they may need updating on
C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      DO 180 I = 1,NDE
 180    SAVE2(I) = SAVE1(I)/YWT(I)
      ETEST = SCNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(REAL(NDE)))
C
C                           The error test failed.  NFAIL keeps track of
C                           multiple failures.  Restore T and the YH
C                           array to their previous values, and prepare
C                           to try the step again.  Compute the optimum
C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL CDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.E0/(BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
          IF (NQ .GT. 1) THEN
            DO 190 I = 1,NDE
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ERDN = SCNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
            RH1 = 1.E0/MAX(1.E0, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
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
          IF (RH*H .EQ. 0.E0) GO TO 400
          CALL CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
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
        CALL CDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (REAL(NSTEP-1)*AVGH + H)/REAL(NSTEP)
      AVGORD = (REAL(NSTEP-1)*AVGORD + REAL(NQ))/REAL(NSTEP)
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
C
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.E0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.E0/REAL(NQ+1)))
            IF (HS .GT. 1.2E0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.E0
              RMAX = RMNORM
              TREND = 1.E0
              CALL CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.E0
              CONVRG = .FALSE.
              CALL CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
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
        RC = 0.E0
        RMAX = RMNORM
        TREND = 1.E0
        CALL CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C                           Consider changing H if NWAIT = 1.  Otherwise
C                           decrease NWAIT by 1.  If NWAIT is then 1 and
C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
C                           a possible order increase on the next step.
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (RH.GT.TRSHLD) CALL CDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
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
          RH1 = 0.E0
        ELSE
          DO 270 I = 1,NDE
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ERDN = SCNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
          RH1 = 1.E0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
        END IF
        RH2 = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.E0
        ELSE
          DO 290 I = 1,NDE
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ERUP = SCNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(REAL(NDE)))
          RH3 = 1.E0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.E0/REAL(NQ+2)))
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
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/REAL(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.E0) RH = MIN(RH, 1.E0/(2.E0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL CDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
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
      END
      SUBROUTINE CDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
C***BEGIN PROLOGUE  CDZRO
C***REFER TO  CDRIV3
C     This is a special purpose version of ZEROIN, modified for use with
C     the CDRIV1 package.
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
C                real result.  This name must be in an
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
C***ROUTINES CALLED  CDNTP
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  CDZRO
      EXTERNAL F
      COMPLEX Y(*), YH(N,*)
      REAL A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND
C***FIRST EXECUTABLE STATEMENT  CDZRO
      ER = 4.E0*UROUND
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
      CMB = 0.5E0*(C - B)
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
      IF (P .LT. 0.E0) THEN
        P = -P
        Q = -Q
      END IF
C                          Update A and check for satisfactory reduction
C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.E0*ACMB .GE. ACBS) THEN
C                                                                 Bisect
          B = 0.5E0*(C + B)
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
        B = 0.5E0*(C + B)
      END IF
C                                             Have completed computation
C                                             for new iterate B.
 20   CALL CDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (FB .EQ. 0.E0) RETURN
      KOUNT = KOUNT + 1
C
C             Decide whether next step is interpolation or extrapolation
C
      IF (SIGN(1.0E0, FB) .EQ. SIGN(1.0E0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A7
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,TRIAD,VECTOR
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
      SUBROUTINE DDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
     8   NDE,NQ,T,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D)
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
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDCOR
      EXTERNAL FA
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
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = DNRM2(N, SAVE2, 1)/SQRT(DBLE(N))
      END IF
      END
      SUBROUTINE DDCST (MAXORD,MINT,ISWFLG,EL,TQ)
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
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDCST
      DOUBLE PRECISION EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
C***FIRST EXECUTABLE STATEMENT  DDCST
      FACTRL(1) = 1.D0
      IF (MAXORD .GE. 2) THEN
        DO 10 I = 2,MAXORD
 10       FACTRL(I) = DBLE(I)*FACTRL(I-1)
      END IF
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
        IF (MAXORD .GE. 3) THEN
          DO 60 J = 3,MAXORD
            EL(2,J) = FACTRL(J-1)
            DO 50 I = 3,J
 50           EL(I,J) = DBLE(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60         EL(J+1,J) = 1.D0
        END IF
        IF (MAXORD .GE. 2) THEN
          DO 80 J = 2,MAXORD
            EL(1,J) = EL(1,J-1) + GAMMA(J)
            EL(2,J) = 1.D0
            DO 80 I = 3,J+1
 80           EL(I,J) = EL(I,J)/(DBLE(I-1)*FACTRL(J-1))
        END IF
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.D0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.D0/GAMMA(J+1)
 100      TQ(3,J) = -1.D0/GAMMA(J+2)
C                                              COMPUTE GEAR COEFFICIENTS
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.D0
        EL(2,1) = 1.D0
        IF (MAXORD .GE. 2) THEN
          DO 130 J = 2,MAXORD
            EL(1,J) = FACTRL(J)
            DO 120 I = 2,J
 120          EL(I,J) = DBLE(J)*EL(I,J-1) + EL(I-1,J-1)
 130        EL(J+1,J) = 1.D0
          SUM = 1.D0
          DO 150 J = 2,MAXORD
            SUM = SUM + 1.D0/DBLE(J)
            DO 150 I = 1,J+1
 150          EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        END IF
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
        IF (MXRD .GE. 2) THEN
          SUM = 1.D0
          DO 200 I = 2,MXRD
            SUM = SUM + 1.D0/DBLE(I)
 200        EL(1+I,1) = -DBLE(I+1)*SUM*GAMMA(I+1)
        END IF
      END IF
      END
      SUBROUTINE DDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,Y,YWT,H,MNTOLD,MTROLD,NFE,RC,YH,
     8   A,CONVRG,EL,IER,IPVT,NQ,NWAIT,RH,RMAX,SAVE2,TQ,TREND,ISWFLG)
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
C***ROUTINES CALLED  DDCST,SDSCL,DGEFA,DGESL,DGBFA,DGBSL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDNTL
      EXTERNAL F, FA
      DOUBLE PRECISION A(MATDIM,*), EL(13,12), EPS, H, HMAX, HNEW, HOLD,
     8     OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX, SMIN,
     8     DNRM2, SUM, SUM0, T, TQ(3,12), TREND, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.D0)
C***FIRST EXECUTABLE STATEMENT  DDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) CALL DDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
        RC = 0.D0
        CONVRG = .FALSE.
        TREND = 1.D0
        RMAX = RMINIT
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              CALL DGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              CALL DGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL DGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            DO 150 I = 1,NDE
              IF(A(I,1) .EQ. 0.D0) THEN
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
          HNEW = H
          RH = H/HOLD
          H = HOLD
          CALL DDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          H = SIGN(MIN(ABS(H), ABS(HNEW)), H)
        END IF
      END IF
      END
      SUBROUTINE DDNTP (H,K,N,NQ,T,TOUT,YH,Y)
C***BEGIN PROLOGUE  DDNTP
C***REFER TO  DDRIV3
C   Subroutine DDNTP interpolates the K-th derivative of Y at TOUT,
C   using the data in the YH array.  If K has a value greater than NQ,
C   the NQ-th derivative is calculated.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDNTP
      DOUBLE PRECISION FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  DDNTP
      KUSED = MIN(K, NQ)
      IF (KUSED .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        FACTOR = 1.D0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*DBLE(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        IF (KUSED .NE. NQ) THEN
          R = ((TOUT - T)/H)
          DO 80 JJ = KUSED+1,NQ
            J = K + 1 + NQ - JJ
            FACTOR = 1.D0
            DO 60 KK = 1,KUSED
 60           FACTOR = FACTOR*DBLE(J-KK)
            DO 70 I = 1,N
 70           Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80         CONTINUE
        END IF
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A4
C***KEYWORDS  BLAS,DOUBLE PRECISION,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
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
     8   NQ,SAVE2,T,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,IER,IPVT,SAVE1,
     8   ISWFLG,BND)
C***BEGIN PROLOGUE  DDPST
C***REFER TO  DDRIV3
C  Subroutine DDPST is called to reevaluate the partials.
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C***ROUTINES CALLED  DGEFA,DGBFA,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDPST
      EXTERNAL F,FA,JACOBN
      DOUBLE PRECISION A(MATDIM,*), BND, DFDY(MATDIM,*), DFDYMX, DY,
     8     EL(13,12), EPSJ, ETA, ETATST, FACTOR, H, SAVE1(*), SAVE2(*),
     8     DNRM2, T, UROUND, Y(*), YH(N,*), YJ, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER, LOOP
      PARAMETER(ETATST = .5D0, ITERMX = 3)
C***FIRST EXECUTABLE STATEMENT  DDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          EPSJ = UROUND**(1.D0/3.D0)
          DO 170 J = 1,N
            DY = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
            IF (DY .EQ. 0.D0) DY = EPSJ
            ITER = 0
 120        YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            Y(J) = YJ
            NFE = NFE + 1
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              DO 130 I = 1,N
                IF (SAVE1(I) .NE. SAVE2(I)) THEN
                  ETA = ABS(SAVE2(I))*UROUND/
     8            (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                  IF (ETA .GE. ETATST) THEN
                    DY = DY*10.D0
                    GO TO 120
                  END IF
                END IF
 130            CONTINUE
            END IF
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = DNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.D0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL DGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          EPSJ = UROUND**(1.D0/3.D0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 265 K = J,N,MW
              DY = EPSJ*MAX(ABS(YWT(K)), ABS(Y(K)))
              IF (DY .EQ. 0.D0) DY = EPSJ
              DFDY(MW,K) = Y(K)
 265          Y(K) = Y(K) + DY
            ITER = 0
 270        CALL F (N, T, Y, SAVE1)
            NFE = NFE + 1
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              LOOP = .FALSE.
              DO 290 K = J,N,MW
                I1 = MAX(1, K-MU)
                I2 = MIN(K+ML, N)
                DO 280 I = I1,I2
                  IF (SAVE1(I) .NE. SAVE2(I)) THEN
                    ETA = ABS(SAVE2(I))*UROUND/
     8              (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                    IF (ETA .GE. ETATST) THEN
                      DY = (Y(K) - DFDY(MW,K))*10.D0
                      Y(K) = DFDY(MW,K) + DY
                      LOOP = .TRUE.
                      GO TO 290
                    END IF
                  END IF
 280              CONTINUE
 290            CONTINUE
              IF (LOOP) GO TO 270
            END IF
            DO 330 K = J,N,MW
              DY = Y(K) - DFDY(MW,K)
              Y(K) = DFDY(MW,K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
 330          CONTINUE
 340        CONTINUE
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
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL DGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
      END IF
      END
      SUBROUTINE DDRIV1 (N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)
C***BEGIN PROLOGUE  DDRIV1
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
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
C           real, independent variable.
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
C             in the call sequence of the user-provided routine F.
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
C             The length of WORK should be at least N*N + 10*N + 225
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to DDRIV1.
C
C***LONG DESCRIPTION
C
C  IV.  USAGE  .........................................................
C
C                   PROGRAM SAMPLE
C                   DOUBLE PRECISION Y(...), WORK(...)
C                   OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C                   N = ...                       Number of equations
C                   T = ...                       Initial point
C                   DO 10 I = 1,N
C              10     Y(I) = ...                  Set initial conditions
C                   TOUT = T
C                   MSTATE = 1
C                   EPS = ...
C                   LENW = ...
C              20   CALL DDRIV1 (N, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C                   IF (MSTATE .GT. 2) STOP
C                   WRITE(6, 100) TOUT, (Y(I), I=1,N)
C                   TOUT = TOUT + 1.
C                   IF (TOUT .LE. 10.) GO TO 20
C              100  FORMAT(...)
C                   END (Sample)
C
C    The user must write a subroutine called F to evaluate the right
C    hand side of the differential equations.  It is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   DOUBLE PRECISION Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C    This computes YDOT = F(Y,T), the right hand side of the
C    differential equations.  Here Y is a vector of length at least N.
C    The actual length of Y is determined by the user's declaration in
C    the program which calls DDRIV1.  Thus the dimensioning of Y in F,
C    while required by FORTRAN convention, does not actually allocate
C    any storage.  When this subroutine is called, the first N
C    components of Y are intermediate approximations to the solution
C    components.  The user should not alter these values.  Here YDOT is
C    a vector of length N.  The user should only compute YDOT(I) for I
C    from 1 to N.
C
C  V.  OTHER COMMUNICATION TO THE USER  ................................
C
C    The solver communicates to the user through the parameters above.
C    In addition it writes diagnostic messages through the standard
C    error handling program XERROR.  That program will terminate the
C    user's run if it detects a probable problem setup error, e.g.,
C    insufficient storage allocated by the user for the WORK array.  For
C    further information see section III-A of the writeup for DDRIV3.
C
C  VI.  REMARKS  .......................................................
C
C    A. There are user-written routines which are only required by
C       DDRIV2 or DDRIV3 when certain parameters are set.  Thus a
C       message warning of unsatisfied externals may be issued during
C       the load or link phase.  This message can be ignored unless it
C       refers to F.
C
C    For other information, see section IV of the writeup for DDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,D1MACH,XERROR
C***END PROLOGUE  DDRIV1
      EXTERNAL F, JACOBN, FA, G
      DOUBLE PRECISION EPS, EWT, G, HMAX, D1MACH, T, TOUT,
     8     WORK(*), Y(*)
      PARAMETER(MXN = 200, IDLIW = 21, MXLIW = IDLIW + MXN)
      INTEGER IWORK(MXLIW)
      CHARACTER MSG*103
      PARAMETER(NROOT = 0, EWT = 1.D0, IERROR = 2, MINT = 2, MITER = 2,
     8          IMPL = 0, ML = 0, MU = 0, MXORD = 5, NDE = 0,
     8          MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  DDRIV1
      IF (N .GT. MXN) THEN
        WRITE(MSG, '(''DDRIV115FE Illegal input.  The number of '',
     8  ''equations,'', I8, '', is greater than the maximum allowed.'')
     8  ') N
        CALL XERROR(MSG, 97, 15, 2)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = SQRT(D1MACH(2))
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 9*N + 204)) THEN
        LNWCHK = N*N + 9*N + 204 + LENIW
        WRITE(MSG, '(''DDRIV116FE Insufficient storage allocated for '',
     8  ''the work array.  The required storage is at least'', I8)')
     8  LNWCHK
        CALL XERROR(MSG, 103, 16, 2)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
          II = I + LENWCM
 20       IWORK(I) = INT(WORK(II))
      END IF
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWT,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      DO 40 I = 1,LENIW
        II = LENWCM + I
 40     WORK(II) = DBLE(IWORK(I))
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE DDRIV2 (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     8   LENW,IWORK,LENIW,G)
C***BEGIN PROLOGUE  DDRIV2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
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
C             G.
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
C             for I from 1 to N.
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
C               N*N + 9*N + 2*NROOT + 204    if MINT is 2, or
C               N*N + 16*N + 2*NROOT + 204   if MINT is 3,
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
C             allocate any storage.
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
C         NSTEP    The number of steps taken.
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
C    B. There are user-written routines which are only required by
C       DDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if it refers to G and NROOT is 0.  A reference to any
C       other unsatisfied external can be ignored.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV2.
C
C  V.  USAGE  ..........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               DOUBLE PRECISION WORK(...), Y(...)           See II. for
C               INTEGER IWORK(...)               required dimensions for
C                                                WORK and IWORK
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C               N = ...                          Number of equations
C               T = ...                          Initial point
C               DO 10 I = 1,N
C          10     Y(I) = ...                     Set initial conditions
C               TOUT = T
C               MSTATE = 1
C               NROOT = 0
C               EPS = ...
C               EWT = ...
C               MINT = 1
C               LENW = ...
C               LENIW = ...
C          20   CALL DDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, G)
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDRIV3,D1MACH,XERROR
C***END PROLOGUE  DDRIV2
      EXTERNAL F, JACOBN, FA, G
      DOUBLE PRECISION EPS, EWT, EWTCOM(1), G, HMAX, D1MACH, T, TOUT,
     8     WORK(*), Y(*)
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, ML = 0, MU = 0, NDE = 0, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  DDRIV2
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(MSG, '(''DDRIV21FE Illegal input.  Improper value for '',
     8  ''the integration method flag,'', I8)') MINT
        CALL XERROR(MSG, 81, 21, 2)
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
      HMAX = SQRT(D1MACH(2))
      CALL DDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE DDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G)
C***BEGIN PROLOGUE  DDRIV3
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
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
C             JACOBN, FA, USERS, and G.
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
C             for I from 1 to N.
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
C                          by the user-written routine named USERS.
C                          This option allows all matrix algebra and
C                          storage decisions to be made by the user.
C                          The routine USERS is called by DDRIV3 when
C                          certain linear systems must be solved.  The
C                          user may choose any method to form, store and
C                          solve these systems in order to obtain the
C                          solution result that is returned to DDRIV3.
C                          In particular, this allows sparse matrix
C                          methods to be used.
C                          The call sequence for this routine is
C
C                           SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2,
C                          8              T, H, EL, IMPL, N, NDE, IFLAG)
C                           DOUBLE PRECISION Y(*), YH(*), YWT(*),
C                          8        SAVE1(*), SAVE2(*), T, H, EL
C
C                          The input variable IFLAG indicates what
C                          action is to be taken.  Subroutine USERS
C                          should perform the following operations,
C                          depending on the value of IFLAG and IMPL.
C
C                          IFLAG = 0
C                            IMPL = 0.  USERS is not called.
C                            IMPL = 1 or 2.  Solve the system
C                                A*X = SAVE2,
C                              returning the result in SAVE2.  The array
C                              SAVE1 can be used as a work array.
C
C                          IFLAG = 1
C                            IMPL = 0.  Compute, decompose and store the
C                              matrix (I - H*EL*J), where I is the
C                              identity matrix and J is the Jacobian
C                              matrix of the right hand side.  The array
C                              SAVE1 can be used as a work array.
C                            IMPL = 1 or 2. Compute, decompose and store
C                              the matrix (A - H*EL*J).  The array SAVE1
C                              can be used as a work array.
C
C                          IFLAG = 2
C                            IMPL = 0.   Solve the system
C                                (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                              returning the result in SAVE2.
C                            IMPL = 1, or 2.  Solve the system
C                              (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                              returning the result in SAVE2.
C                            The array SAVE1 should not be altered.
C
C                          When using a value of MITER = 3, the
C                          subroutine FA is not required, even if IMPL
C                          is not 0.  For further information on using
C                          this option, see section IV-F below.
C
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
C         1,2  N*N+(MXORD+4)*N     2*N*N+(MXORD+4)*N   N*N+(MXORD+5)*N
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204
C
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +
C              (MXORD+5)*N +       (MXORD+6)*N +       (MXORD+6)*N +
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
C             K=I-J+MU+1.
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
C             and Y.
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
C             allocate any storage.
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
C         9   Fatal        Improper value for MINT, MITER and/or IMPL.
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
C         NSTEP    The number of steps taken.
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
C    C. There are user-written routines which are only required by
C       DDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if: it refers to JACOBN and MITER is not 1 or 4, or
C       it refers to FA and IMPL is 0 or MITER is 3, or it refers to
C       USERS and MITER is not 3, or it refers to G and NROOT is 0.
C
C    D. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to DDRIV3.
C
C    E. Changing parameters during an integration.
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
C    F. As the price for complete control of matrix algebra, the DDRIV3
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
C              EPSJ = UROUND**(1.D0/3.D0)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.D0) R = EPSJ
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
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
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  DDSTP,SDNTP,SDZRO,DGEFA,DGESL,DGBFA,DGBSL,DNRM2,
C                    D1MACH,XERROR
C***END PROLOGUE  DDRIV3
      EXTERNAL F, JACOBN, FA, G
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
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20)
      PARAMETER(INDPRT = 21, INDPVT = 22)
C***FIRST EXECUTABLE STATEMENT  DDRIV3
      UROUND = D1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = D1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.D0) THEN
        WRITE(MSG, '(''DDRIV36FE Illegal input.  EPS,'', D16.8,
     8  '', is negative.'')') EPS
        CALL XERROR(MSG, 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''DDRIV37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
        CALL XERROR(MSG, 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''DDRIV314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
        CALL XERROR(MSG, 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1)) THEN
        WRITE(MSG, '(''DDRIV39FE Illegal input.  Improper value for '',
     8  ''MINT, MITER and/or IMPL.'')')
        CALL XERROR(MSG, 69, 9, 2)
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
        CALL XERROR(MSG, 164, 10, 2)
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
      IA = ITROOT + NROOT
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
        CALL XERROR(MSG, 162, 8, 2)
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
        IF (T .EQ. TOUT) RETURN
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
        WORK(IAVGRD) = 0.D0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
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
      IF(NROOT .NE. 0) THEN
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
          CALL XERROR(MSG, 124, 2, 0)
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
          CALL XERROR(MSG, 124, 2, 0)
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
 200      WORK(JGNOW) = G (N, T, Y, I)
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
          CALL F (N, T, Y, WORK(ISAVE2))
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y,WORK(IYH),WORK(IYWT),WORK(ISAVE1),WORK(ISAVE2),
     8                 T, H, WORK(IEL), IMPL, N, NDECOM, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              CALL DGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (N, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              CALL DGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL DGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF(WORK(JA) .EQ. 0.D0) GO TO 690
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
        CALL XERROR(MSG, 152, 4, 1)
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
        CALL XERROR(MSG, 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''DDRIV33WRN At T,'', D16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', D16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG, 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL DDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL DDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, N,
     8            NDECOM, WORK(IYWT), UROUND, WORK(IAVGH), WORK(IAVGRD),
     8            WORK(IH), WORK(IHUSED), IWORK(IJTASK), IWORK(IMNTLD),
     8            IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IHOLD),
     8            IWORK(INDPVT), JSTATE, IWORK(INQ), IWORK(INWAIT),
     8            WORK(IRC), WORK(IRMAX), WORK(ISAVE1), WORK(ISAVE2),
     8            WORK(ITQ), WORK(ITREND), MINT, IWORK(IMTRSV),
     8            IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (ABS(H) .GE. UROUND*ABS(T) .AND. NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (N, T, Y, I)
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
                TLAST = T - WORK(IHUSED)
                IROOT = I
                TROOT = T
                CALL DDZRO (AE, G, H, N, IWORK(INQ), IROOT, RE, T,
     8                      WORK(IYH), UROUND,  TROOT, TLAST,
     8                      WORK(JGNOW), GLAST,  WORK(ISAVE1))
                WORK(JTROOT) = TROOT
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
C                                        Fatal errors are processed here
C
 670  WRITE(MSG, '(''DDRIV311FE At T,'', D16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL XERROR(MSG, 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''DDRIV312FE At T,'', D16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG, 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''DDRIV313FE At T,'', D16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG, 74, 13, 2)
      RETURN
      END
      SUBROUTINE DDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)
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
     8   MITER,ML,MU,N,NDE,YWT,UROUND,AVGH,AVGORD,H,HUSED,JTASK,MNTOLD,
     8   MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,HOLD,IPVT,
     8   JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,MTRSV,
     8   MXRDSV)
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
C***ROUTINES CALLED  DDNTL,SDPST,SDCOR,SDPSC,SDSCL,DNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860513   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDSTP
      EXTERNAL F, JACOBN, FA
      DOUBLE PRECISION A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, H, HMAX, HN, HOLD, HS, HUSED, NUMER, RC,
     8     RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM, SAVE1(*),
     8     SAVE2(*), DNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, UROUND,
     8     Y(*), YH(N,*), YWT(*), Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3D0, BIAS2 = 1.2D0, BIAS3 = 1.4D0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3D0, RMFAIL = 2.D0,
     8          RMNORM = 10.D0, TRSHLD = 1.D0)
C***FIRST EXECUTABLE STATEMENT  DDSTP
      BND = 0.D0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL DDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
        IF (H .EQ. 0.D0) GO TO 400
        IF (IER) GO TO 420
      END IF
      IER = .FALSE.
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
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL DDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, IER, IPVT, SAVE1, ISWFLG, BND)
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
     8            MU, N, NDE, NQ, T, Y, YH, YWT,  EVALFA, SAVE1,
     8            SAVE2,  A, D)
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
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
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
      END
      SUBROUTINE DDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
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
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  DDZRO
      EXTERNAL F
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
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
      INTEGER LDA,N,ML,MU,IPVT(1),INFO
      DOUBLE PRECISION ABD(LDA,1)
C
C     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.
C
C     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
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
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C     FORTRAN MAX0,MIN0
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C
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
      INTEGER LDA,N,ML,MU,IPVT(1),JOB
      DOUBLE PRECISION ABD(LDA,1),B(1)
C
C     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGBCO OR DGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0
C        OR DGBFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN MIN0
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C
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
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
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
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR DGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C
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
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
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
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A6
C***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
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
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
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
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D1A2
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAXIMUM COMPONENT,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find largest component of d.p. vector
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
      DOUBLE PRECISION DX(1),DMAX,XMAG
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
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
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
      SUBROUTINE SDCOR (DFDY,EL,FA,H,IMPL,IPVT,MATDIM,MITER,ML,MU,N,
     8   NDE,NQ,T,Y,YH,YWT,EVALFA,SAVE1,SAVE2,A,D)
C***BEGIN PROLOGUE  SDCOR
C***REFER TO  SDRIV3
C  Subroutine SDCOR is called to compute corrections to the Y array.
C  In the case of functional iteration, update Y directly from the
C  result of the last call to F.
C  In the case of the chord method, compute the corrector error and
C  solve the linear system with that as right hand side and DFDY as
C  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4,
C  or 5.
C***ROUTINES CALLED  SGESL,SGBSL,SNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDCOR
      EXTERNAL FA
      REAL A(MATDIM,*), D, DFDY(MATDIM,*), EL(13,12), H,
     8     SAVE1(*), SAVE2(*), SNRM2, T, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL EVALFA
C***FIRST EXECUTABLE STATEMENT  SDCOR
      IF (MITER .EQ. 0) THEN
        DO 100 I = 1,N
 100      SAVE1(I) = (H*SAVE2(I) - YH(I,2) - SAVE1(I))/YWT(I)
        D = SNRM2(N, SAVE1, 1)/SQRT(REAL(N))
        DO 105 I = 1,N
 105      SAVE1(I) = H*SAVE2(I) - YH(I,2)
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (IMPL .EQ. 0) THEN
          DO 130 I = 1,N
 130        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
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
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 180 I = 1,N
 180        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL SGESL (DFDY, MATDIM, N, IPVT, SAVE2, 0)
        DO 200 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 200      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (IMPL .EQ. 0) THEN
          DO 230 I = 1,N
 230        SAVE2(I) = H*SAVE2(I) - YH(I,2) - SAVE1(I)
        ELSE IF (IMPL .EQ. 1) THEN
          IF (EVALFA) THEN
            CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
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
          ELSE
            EVALFA = .TRUE.
          END IF
          DO 280 I = 1,N
 280        SAVE2(I) = H*SAVE2(I) - A(I,1)*(YH(I,2) + SAVE1(I))
        END IF
        CALL SGBSL (DFDY, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
        DO 300 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 300      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 2
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
        DO 320 I = 1,N
          SAVE1(I) = SAVE1(I) + SAVE2(I)
 320      SAVE2(I) = SAVE2(I)/YWT(I)
        D = SNRM2(N, SAVE2, 1)/SQRT(REAL(N))
      END IF
      END
      SUBROUTINE SDCST (MAXORD,MINT,ISWFLG,EL,TQ)
C***BEGIN PROLOGUE  SDCST
C***REFER TO  SDRIV3
C  SDCST is called by SDNTL and sets coefficients used by the core
C  integrator SDSTP.  The array EL determines the basic method.
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
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDCST
      REAL EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
C***FIRST EXECUTABLE STATEMENT  SDCST
      FACTRL(1) = 1.E0
      IF (MAXORD .GE. 2) THEN
        DO 10 I = 2,MAXORD
 10       FACTRL(I) = REAL(I)*FACTRL(I-1)
      END IF
C                                             COMPUTE ADAMS COEFFICIENTS
      IF (MINT .EQ. 1) THEN
        GAMMA(1) = 1.E0
        DO 40 I = 1,MAXORD+1
          SUM = 0.E0
          DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/REAL(I-J+2)
 40       GAMMA(I+1) = SUM
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        EL(2,2) = 1.E0
        EL(3,2) = 1.E0
        IF (MAXORD .GE. 3) THEN
          DO 60 J = 3,MAXORD
            EL(2,J) = FACTRL(J-1)
            DO 50 I = 3,J
 50           EL(I,J) = REAL(J-1)*EL(I,J-1) + EL(I-1,J-1)
 60         EL(J+1,J) = 1.E0
        END IF
        IF (MAXORD .GE. 2) THEN
          DO 80 J = 2,MAXORD
            EL(1,J) = EL(1,J-1) + GAMMA(J)
            EL(2,J) = 1.E0
            DO 80 I = 3,J+1
 80           EL(I,J) = EL(I,J)/(REAL(I-1)*FACTRL(J-1))
        END IF
        DO 100 J = 1,MAXORD
          TQ(1,J) = -1.E0/(FACTRL(J)*GAMMA(J))
          TQ(2,J) = -1.E0/GAMMA(J+1)
 100      TQ(3,J) = -1.E0/GAMMA(J+2)
C                                              COMPUTE GEAR COEFFICIENTS
      ELSE IF (MINT .EQ. 2) THEN
        EL(1,1) = 1.E0
        EL(2,1) = 1.E0
        IF (MAXORD .GE. 2) THEN
          DO 130 J = 2,MAXORD
            EL(1,J) = FACTRL(J)
            DO 120 I = 2,J
 120          EL(I,J) = REAL(J)*EL(I,J-1) + EL(I-1,J-1)
 130        EL(J+1,J) = 1.E0
          SUM = 1.E0
          DO 150 J = 2,MAXORD
            SUM = SUM + 1.E0/REAL(J)
            DO 150 I = 1,J+1
 150          EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
        END IF
        DO 170 J = 1,MAXORD
          IF (J .GT. 1) TQ(1,J) = 1.E0/FACTRL(J-1)
          TQ(2,J) = REAL(J+1)/EL(1,J)
 170      TQ(3,J) = REAL(J+2)/EL(1,J)
      END IF
C                          Compute constants used in the stiffness test.
C                          These are the ratio of TQ(2,NQ) for the Gear
C                          methods to those for the Adams methods.
      IF (ISWFLG .EQ. 3) THEN
        MXRD = MIN(MAXORD, 5)
        IF (MINT .EQ. 2) THEN
          GAMMA(1) = 1.E0
          DO 190 I = 1,MXRD
            SUM = 0.E0
            DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/REAL(I-J+2)
 190        GAMMA(I+1) = SUM
        END IF
        IF (MXRD .GE. 2) THEN
          SUM = 1.E0
          DO 200 I = 2,MXRD
            SUM = SUM + 1.E0/REAL(I)
 200        EL(1+I,1) = -REAL(I+1)*SUM*GAMMA(I+1)
        END IF
      END IF
      END
      SUBROUTINE SDNTL (EPS,F,FA,HMAX,HOLD,IMPL,JTASK,MATDIM,MAXORD,
     8   MINT,MITER,ML,MU,N,NDE,SAVE1,T,Y,YWT,H,MNTOLD,MTROLD,NFE,RC,YH,
     8   A,CONVRG,EL,IER,IPVT,NQ,NWAIT,RH,RMAX,SAVE2,TQ,TREND,ISWFLG)
C***BEGIN PROLOGUE  SDNTL
C***REFER TO  SDRIV3
C  Subroutine SDNTL is called to set parameters on the first call
C  to SDSTP, on an internal restart, or when the user has altered
C  MINT, MITER, and/or H.
C  On the first call, the order is set to 1 and the initial derivatives
C  are calculated.  RMAX is the maximum ratio by which H can be
C  increased in one step.  It is initially RMINIT to compensate
C  for the small initial H, but then is normally equal to RMNORM.
C  If a failure occurs (in corrector convergence or error test), RMAX
C  is set at RMFAIL for the next increase.
C  If the caller has changed MINT, or if JTASK = 0, SDCST is called
C  to set the coefficients of the method.  If the caller has changed H,
C  YH must be rescaled.  If H or MINT has been changed, NWAIT is
C  reset to NQ + 2 to prevent further increases in H for that many
C  steps.  Also, RC is reset.  RC is the ratio of new to old values of
C  the coefficient L(0)*H.  If the caller has changed MITER, RC is
C  set to 0 to force the partials to be updated, if partials are used.
C***ROUTINES CALLED  SDCST,SDSCL,SGEFA,SGESL,SGBFA,SGBSL,SNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDNTL
      EXTERNAL F,FA
      REAL A(MATDIM,*), EL(13,12), EPS, H, HMAX, HNEW, HOLD,
     8     OLDL0, RC, RH, RMAX, RMINIT, SAVE1(*), SAVE2(*), SMAX, SMIN,
     8     SNRM2, SUM, SUM0, T, TQ(3,12), TREND, Y(*), YH(N,*), YWT(*)
      INTEGER IPVT(*)
      LOGICAL CONVRG, IER
      PARAMETER(RMINIT = 10000.E0)
C***FIRST EXECUTABLE STATEMENT  SDNTL
      IER = .FALSE.
      IF (JTASK .GE. 0) THEN
        IF (JTASK .EQ. 0) CALL SDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
        RC = 0.E0
        CONVRG = .FALSE.
        TREND = 1.E0
        RMAX = RMINIT
        NQ = 1
        NWAIT = 3
        CALL F (N, T, Y, SAVE2)
        NFE = NFE + 1
        IF (IMPL .NE. 0) THEN
          IF (MITER .EQ. 3) THEN
            IFLAG = 0
            CALL USERS (Y, YH, YWT, SAVE1, SAVE2, T, H, EL, IMPL, N,
     8                  NDE, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
              CALL SGEFA (A, MATDIM, N, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL SGESL (A, MATDIM, N, IPVT, SAVE2, 0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
              CALL SGBFA (A, MATDIM, N, ML, MU, IPVT, INFO)
              IF (INFO .NE. 0) THEN
                IER = .TRUE.
                RETURN
              END IF
              CALL SGBSL (A, MATDIM, N, ML, MU, IPVT, SAVE2, 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
            DO 150 I = 1,NDE
              IF(A(I,1) .EQ. 0.E0) THEN
                IER = .TRUE.
                RETURN
              ELSE
                SAVE2(I) = SAVE2(I)/A(I,1)
              END IF
 150          CONTINUE
            DO 155 I = NDE+1,N
 155          A(I,1) = 0.E0
          END IF
        END IF
        DO 170 I = 1,NDE
 170      SAVE1(I) = SAVE2(I)/YWT(I)
        SUM = SNRM2(NDE, SAVE1, 1)
        SUM0 = 1.E0/MAX(1.E0, ABS(T))
        SMAX = MAX(SUM0, SUM)
        SMIN = MIN(SUM0, SUM)
        SUM = SMAX*SQRT(1.E0 + (SMIN/SMAX)**2)/SQRT(REAL(NDE))
        H = SIGN(MIN(2.E0*EPS/SUM, ABS(H)), H)
        DO 180 I = 1,N
 180      YH(I,2) = H*SAVE2(I)
      ELSE
        IF (MITER .NE. MTROLD) THEN
          MTROLD = MITER
          RC = 0.E0
          CONVRG = .FALSE.
        END IF
        IF (MINT .NE. MNTOLD) THEN
          MNTOLD = MINT
          OLDL0 = EL(1,NQ)
          CALL SDCST (MAXORD, MINT, ISWFLG,  EL, TQ)
          RC = RC*EL(1,NQ)/OLDL0
          NWAIT = NQ + 2
        END IF
        IF (H .NE. HOLD) THEN
          NWAIT = NQ + 2
          HNEW = H
          RH = H/HOLD
          H = HOLD
          CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
          H = SIGN(MIN(ABS(H), ABS(HNEW)), H)
        END IF
      END IF
      END
      SUBROUTINE SDNTP (H,K,N,NQ,T,TOUT,YH,Y)
C***BEGIN PROLOGUE  SDNTP
C***REFER TO  SDRIV3
C   Subroutine SDNTP interpolates the K-th derivative of Y at TOUT,
C   using the data in the YH array.  If K has a value greater than NQ,
C   the NQ-th derivative is calculated.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDNTP
      REAL FACTOR, H, R, T, TOUT, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  SDNTP
      KUSED = MIN(K, NQ)
      IF (KUSED .EQ. 0) THEN
        DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
        R = ((TOUT - T)/H)
        DO 20 JJ = 1,NQ
          J = NQ + 1 - JJ
          DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
      ELSE
        FACTOR = 1.E0
        DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*REAL(NQ+1-KK)
        DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
        IF (KUSED .NE. NQ) THEN
          R = ((TOUT - T)/H)
          DO 80 JJ = KUSED+1,NQ
            J = K + 1 + NQ - JJ
            FACTOR = 1.E0
            DO 60 KK = 1,KUSED
 60           FACTOR = FACTOR*REAL(J-KK)
            DO 70 I = 1,N
 70           Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80         CONTINUE
        END IF
        DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
      END IF
      END
      SUBROUTINE SDPSC (KSGN,N,NQ,YH)
C***BEGIN PROLOGUE  SDPSC
C***REFER TO  SDRIV3
C     This subroutine computes the predicted YH values by effectively
C     multiplying the YH array by the Pascal triangle matrix when KSGN
C     is +1, and performs the inverse function when KSGN is -1.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDPSC
      REAL YH(N,*)
C***FIRST EXECUTABLE STATEMENT  SDPSC
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
      SUBROUTINE SDPST (EL,F,FA,H,IMPL,JACOBN,MATDIM,MITER,ML,MU,N,NDE,
     8   NQ,SAVE2,T,Y,YH,YWT,UROUND,NFE,NJE,A,DFDY,IER,IPVT,SAVE1,
     8   ISWFLG,BND)
C***BEGIN PROLOGUE  SDPST
C***REFER TO  SDRIV3
C  Subroutine SDPST is called to reevaluate the partials.
C  If MITER is 1, 2, 4, or 5, the matrix
C  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU
C  decomposition, with the results also stored in DFDY.
C***ROUTINES CALLED  SGEFA,SGBFA,SNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850320   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDPST
      EXTERNAL F,FA,JACOBN
      REAL A(MATDIM,*), BND, DFDY(MATDIM,*), DFDYMX, DY,
     8     EL(13,12), EPSJ, ETA, ETATST, FACTOR, H, SAVE1(*), SAVE2(*),
     8     SNRM2, T, UROUND, Y(*), YH(N,*), YJ, YWT(*)
      INTEGER IPVT(*)
      LOGICAL IER, LOOP
      PARAMETER(ETATST = .5E0, ITERMX = 3)
C***FIRST EXECUTABLE STATEMENT  SDPST
      NJE = NJE + 1
      IER = .FALSE.
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        IF (MITER .EQ. 1) THEN
          CALL JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
          IF (ISWFLG .EQ. 3) BND = SNRM2(N*N, DFDY, 1)
          FACTOR = -EL(1,NQ)*H
          DO 110 J = 1,N
            DO 110 I = 1,N
 110          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 2) THEN
          EPSJ = UROUND**(1.E0/3.E0)
          DO 170 J = 1,N
            DY = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
            IF (DY .EQ. 0.E0) DY = EPSJ
            ITER = 0
 120        YJ = Y(J)
            Y(J) = Y(J) + DY
            CALL F (N, T, Y, SAVE1)
            Y(J) = YJ
            NFE = NFE + 1
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              DO 130 I = 1,N
                IF (SAVE1(I) .NE. SAVE2(I)) THEN
                  ETA = ABS(SAVE2(I))*UROUND/
     8            (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                  IF (ETA .GE. ETATST) THEN
                    DY = DY*10.E0
                    GO TO 120
                  END IF
                END IF
 130            CONTINUE
            END IF
            FACTOR = -EL(1,NQ)*H/DY
            DO 140 I = 1,N
 140          DFDY(I,J) = (SAVE1(I) - SAVE2(I))*FACTOR
 170        CONTINUE
          IF (ISWFLG .EQ. 3) BND = SNRM2(N*N, DFDY, 1)/(-EL(1,NQ)*H)
        END IF
        IF (IMPL .EQ. 0) THEN
          DO 190 I = 1,N
 190        DFDY(I,I) = DFDY(I,I) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 210 J = 1,N
            DO 210 I = 1,N
 210          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 230 I = 1,NDE
 230        DFDY(I,I) = DFDY(I,I) + A(I,1)
        END IF
        CALL SGEFA (DFDY, MATDIM, N, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        IF (MITER .EQ. 4) THEN
          CALL JACOBN (N, T, Y, DFDY(ML+1,1), MATDIM, ML, MU)
          FACTOR = -EL(1,NQ)*H
          MW = ML + MU + 1
          DO 260 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 260 I = I1,I2
 260          DFDY(I,J) = FACTOR*DFDY(I,J)
        ELSE IF (MITER .EQ. 5) THEN
          EPSJ = UROUND**(1.E0/3.E0)
          MW = ML + MU + 1
          J2 = MIN(MW, N)
          DO 340 J = 1,J2
            DO 265 K = J,N,MW
              DY = EPSJ*MAX(ABS(YWT(K)), ABS(Y(K)))
              IF (DY .EQ. 0.E0) DY = EPSJ
              DFDY(MW,K) = Y(K)
 265          Y(K) = Y(K) + DY
            ITER = 0
 270        CALL F (N, T, Y, SAVE1)
            NFE = NFE + 1
            ITER = ITER + 1
            IF (ITER .LT. ITERMX) THEN
              LOOP = .FALSE.
              DO 290 K = J,N,MW
                I1 = MAX(1, K-MU)
                I2 = MIN(K+ML, N)
                DO 280 I = I1,I2
                  IF (SAVE1(I) .NE. SAVE2(I)) THEN
                    ETA = ABS(SAVE2(I))*UROUND/
     8              (ABS(SAVE2(I) - SAVE1(I)) + ABS(SAVE2(I))*UROUND)
                    IF (ETA .GE. ETATST) THEN
                      DY = (Y(K) - DFDY(MW,K))*10.E0
                      Y(K) = DFDY(MW,K) + DY
                      LOOP = .TRUE.
                      GO TO 290
                    END IF
                  END IF
 280              CONTINUE
 290            CONTINUE
              IF (LOOP) GO TO 270
            END IF
            DO 330 K = J,N,MW
              DY = Y(K) - DFDY(MW,K)
              Y(K) = DFDY(MW,K)
              FACTOR = -EL(1,NQ)*H/DY
              I1 = MAX(ML+1, MW+1-K)
              I2 = MIN(MW+N-K, MW+ML)
              DO 300 I = I1,I2
                I3 = K + I - MW
 300            DFDY(I,K) = FACTOR*(SAVE1(I3) - SAVE2(I3))
 330          CONTINUE
 340        CONTINUE
        END IF
        IF (ISWFLG .EQ. 3) THEN
          DFDYMX = 0.E0
          DO 345 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 345 I = I1,I2
 345          DFDYMX = MAX(DFDYMX, ABS(DFDY(I,J)))
          BND = 0.E0
          IF (DFDYMX .NE. 0.E0) THEN
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
 360        DFDY(MW,J) = DFDY(MW,J) + 1.E0
        ELSE IF (IMPL .EQ. 1) THEN
          CALL FA (N, T, Y, A(ML+1,1), MATDIM, ML, MU, NDE)
          DO 380 J = 1,N
            I1 = MAX(ML+1, MW+1-J)
            I2 = MIN(MW+N-J, MW+ML)
            DO 380 I = I1,I2
 380          DFDY(I,J) = DFDY(I,J) + A(I,J)
        ELSE IF (IMPL .EQ. 2) THEN
          CALL FA (N, T, Y, A, MATDIM, ML, MU, NDE)
          DO 400 J = 1,NDE
 400        DFDY(MW,J) =  DFDY(MW,J) + A(J,1)
        END IF
        CALL SGBFA (DFDY, MATDIM, N, ML, MU, IPVT, INFO)
        IF (INFO .NE. 0) IER = .TRUE.
      ELSE IF (MITER .EQ. 3) THEN
        IFLAG = 1
        CALL USERS (Y, YH(1,2), YWT, SAVE1, SAVE2, T, H, EL(1,NQ), IMPL,
     8              N, NDE, IFLAG)
      END IF
      END
      SUBROUTINE SDRIV1 (N,T,Y,TOUT,MSTATE,EPS,WORK,LENW)
C***BEGIN PROLOGUE  SDRIV1
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             SINGLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of SDRIV1 is to solve N (200 or fewer)
C            ordinary differential equations of the form
C            dY(I)/dT = F(Y(I),T), given the initial conditions
C            Y(I) = YI.  SDRIV1 uses single precision arithmetic.
C***DESCRIPTION
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
C           real, independent variable.
C
C    As an aid in selecting the proper program, the following is a
C    discussion of the important options or restrictions associated with
C    each program:
C
C      A. SDRIV1 should be tried first for those routine problems with
C         no more than 200 differential equations.  Internally this
C         routine has two important technical defaults:
C           1. Numerical approximation of the Jacobian matrix of the
C              right hand side is used.
C           2. The stiff solver option is used.
C         Most users of SDRIV1 should not have to concern themselves
C         with these details.
C
C      B. SDRIV2 should be considered for those problems for which
C         SDRIV1 is inadequate (SDRIV2 has no explicit restriction on
C         the number of differential equations.)  For example, SDRIV1
C         may have difficulty with problems having zero initial
C         conditions and zero derivatives.  In this case SDRIV2, with an
C         appropriate value of the parameter EWT, should perform more
C         efficiently.  SDRIV2 provides three important additional
C         options:
C           1. The nonstiff equation solver (as well as the stiff
C              solver) is available.
C           2. The root-finding option is available.
C           3. The program can dynamically select either the non-stiff
C              or the stiff methods.
C         Internally this routine also defaults to the numerical
C         approximation of the Jacobian matrix of the right hand side.
C
C      C. SDRIV3 is the most flexible, and hence the most complex, of
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
C    The function of SDRIV1 is to solve N (200 or fewer) ordinary
C    differential equations of the form dY(I)/dT = F(Y(I),T), given the
C    initial conditions Y(I) = YI.  SDRIV1 is to be called once for each
C    output point.
C
C  III.  PARAMETERS  ...................................................
C
C    The user should use parameter names in the call sequence of SDRIV1
C    for those quantities whose value may be altered by SDRIV1.  The
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
C             in the call sequence of the user-provided routine F.
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
C                  user.  Unless SDRIV1 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, SDRIV1 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling SDRIV1
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling SDRIV1
C                  again.
C
C    EPS    = On input, the requested relative accuracy in all solution
C             components.  On output, the adjusted relative accuracy if
C             the input value was too small.  The value of EPS should be
C             set as large as is reasonable, because the amount of work
C             done by SDRIV1 increases as EPS decreases.
C
C    WORK
C    LENW   = (Input)
C             WORK is an array of LENW real words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       REAL WORK(...)
C             The length of WORK should be at least N*N + 10*N + 225
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to SDRIV1.
C
C***LONG DESCRIPTION
C
C  IV.  USAGE  .........................................................
C
C                   PROGRAM SAMPLE
C                   REAL Y(...), WORK(...)
C                   OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C                   N = ...                       Number of equations
C                   T = ...                       Initial point
C                   DO 10 I = 1,N
C              10     Y(I) = ...                  Set initial conditions
C                   TOUT = T
C                   MSTATE = 1
C                   EPS = ...
C                   LENW = ...
C              20   CALL SDRIV1 (N, T, Y, TOUT, MSTATE, EPS, WORK, LENW)
C                   IF (MSTATE .GT. 2) STOP
C                   WRITE(6, 100) TOUT, (Y(I), I=1,N)
C                   TOUT = TOUT + 1.
C                   IF (TOUT .LE. 10.) GO TO 20
C              100  FORMAT(...)
C                   END (Sample)
C
C    The user must write a subroutine called F to evaluate the right
C    hand side of the differential equations.  It is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   REAL Y(*), YDOT(*)
C                     .
C                     .
C                   YDOT(1) = ...
C                     .
C                     .
C                   YDOT(N) = ...
C                   END (Sample)
C    This computes YDOT = F(Y,T), the right hand side of the
C    differential equations.  Here Y is a vector of length at least N.
C    The actual length of Y is determined by the user's declaration in
C    the program which calls SDRIV1.  Thus the dimensioning of Y in F,
C    while required by FORTRAN convention, does not actually allocate
C    any storage.  When this subroutine is called, the first N
C    components of Y are intermediate approximations to the solution
C    components.  The user should not alter these values.  Here YDOT is
C    a vector of length N.  The user should only compute YDOT(I) for I
C    from 1 to N.
C
C  V.  OTHER COMMUNICATION TO THE USER  ................................
C
C    The solver communicates to the user through the parameters above.
C    In addition it writes diagnostic messages through the standard
C    error handling program XERROR.  That program will terminate the
C    user's run if it detects a probable problem setup error, e.g.,
C    insufficient storage allocated by the user for the WORK array.  For
C    further information see section III-A of the writeup for SDRIV3.
C
C  VI.  REMARKS  .......................................................
C
C    A. There are user-written routines which are only required by
C       SDRIV2 or SDRIV3 when certain parameters are set.  Thus a
C       message warning of unsatisfied externals may be issued during
C       the load or link phase.  This message can be ignored unless it
C       refers to F.
C
C    For other information, see section IV of the writeup for SDRIV3.
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  SDRIV3,R1MACH,XERROR
C***END PROLOGUE  SDRIV1
      EXTERNAL F, JACOBN, FA, G
      REAL EPS, EWT, G, HMAX, R1MACH, T, TOUT,
     8     WORK(*), Y(*)
      PARAMETER(MXN = 200, IDLIW = 21, MXLIW = IDLIW + MXN)
      INTEGER IWORK(MXLIW)
      CHARACTER MSG*103
      PARAMETER(NROOT = 0, EWT = 1.E0, IERROR = 2, MINT = 2, MITER = 2,
     8          IMPL = 0, ML = 0, MU = 0, MXORD = 5, NDE = 0,
     8          MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  SDRIV1
      IF (N .GT. MXN) THEN
        WRITE(MSG, '(''SDRIV115FE Illegal input.  The number of '',
     8  ''equations,'', I8, '', is greater than the maximum allowed.'')
     8  ') N
        CALL XERROR(MSG, 97, 15, 2)
        RETURN
      END IF
      IF (MSTATE .GT. 0) THEN
        NSTATE = MSTATE
        NTASK = 1
      ELSE
        NSTATE = - MSTATE
        NTASK = 3
      END IF
      HMAX = SQRT(R1MACH(2))
      LENIW = N + IDLIW
      LENWCM = LENW - LENIW
      IF (LENWCM .LT. (N*N + 9*N + 204)) THEN
        LNWCHK = N*N + 9*N + 204 + LENIW
        WRITE(MSG, '(''SDRIV116FE Insufficient storage allocated for '',
     8  ''the work array.  The required storage is at least'', I8)')
     8  LNWCHK
        CALL XERROR(MSG, 103, 16, 2)
        RETURN
      END IF
      IF (NSTATE .NE. 1) THEN
        DO 20 I = 1,LENIW
          II = I + LENWCM
 20       IWORK(I) = INT(WORK(II))
      END IF
      CALL SDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWT,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENWCM, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      DO 40 I = 1,LENIW
        II = LENWCM + I
 40     WORK(II) = REAL(IWORK(I))
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE SDRIV2 (N,T,Y,F,TOUT,MSTATE,NROOT,EPS,EWT,MINT,WORK,
     8   LENW,IWORK,LENIW,G)
C***BEGIN PROLOGUE  SDRIV2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             SINGLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of SDRIV2 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  SDRIV2 uses single precision arithmetic.
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The function of SDRIV2 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.
C    SDRIV2 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of SDRIV2
C    for those quantities whose value may be altered by SDRIV2.  The
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
C             G.
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   REAL Y(*), YDOT(*)
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
C             user's declaration in the program which calls SDRIV2.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.
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
C                  user.  Unless SDRIV2 is to be reinitialized, only the
C                  sign of MSTATE may be changed by the user.  (As a
C                  convenience to the user who may wish to put out the
C                  initial conditions, SDRIV2 can be called with
C                  MSTATE=+1(-1), and TOUT=T.  In this case the program
C                  will return with MSTATE unchanged, i.e.,
C                  MSTATE=+1(-1).)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  1000 steps without reaching TOUT.  The user can
C                  continue the integration by simply calling SDRIV2
C                  again.  Other than an error in problem setup, the
C                  most likely cause for this condition is trying to
C                  integrate a stiff set of equations with the non-stiff
C                  integrator option. (See description of MINT below.)
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling SDRIV2
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling SDRIV2 again.
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  SDRIV2 attempts
C             to find the value of T at which one of the equations
C             changes sign.  SDRIV2 can find at most one root per
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
C             because the amount of work done by SDRIV2 increases as
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
C             WORK is an array of LENW real words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       REAL WORK(...)
C             The length of WORK should be at least
C               16*N + 2*NROOT + 204         if MINT is 1, or
C               N*N + 9*N + 2*NROOT + 204    if MINT is 2, or
C               N*N + 16*N + 2*NROOT + 204   if MINT is 3,
C             and LENW should be set to the value used.  The contents of
C             WORK should not be disturbed between calls to SDRIV2.
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
C             of IWORK should not be disturbed between calls to SDRIV2.
C
C    G      = A real FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to
C             obtain the value of each of the NROOT equations for which
C             a root is desired.  G is of the form:
C                   REAL FUNCTION G (N, T, Y, IROOT)
C                   REAL Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls SDRIV2.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.
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
C         NSTEP    The number of steps taken.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. On any return from SDRIV2 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    B. There are user-written routines which are only required by
C       SDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if it refers to G and NROOT is 0.  A reference to any
C       other unsatisfied external can be ignored.
C
C    C. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to SDRIV2.
C
C  V.  USAGE  ..........................................................
C
C               PROGRAM SAMPLE
C               EXTERNAL F
C               REAL WORK(...), Y(...)           See II. for
C               INTEGER IWORK(...)               required dimensions for
C                                                WORK and IWORK
C               OPEN(FILE='TAPE6', UNIT=6, STATUS='NEW')
C               N = ...                          Number of equations
C               T = ...                          Initial point
C               DO 10 I = 1,N
C          10     Y(I) = ...                     Set initial conditions
C               TOUT = T
C               MSTATE = 1
C               NROOT = 0
C               EPS = ...
C               EWT = ...
C               MINT = 1
C               LENW = ...
C               LENIW = ...
C          20   CALL SDRIV2 (N, T, Y, F, TOUT, MSTATE, NROOT, EPS, EWT,
C              8             MINT, WORK, LENW, IWORK, LENIW, G)
C               IF (MSTATE .GT. 2) STOP
C               WRITE(6, 100) TOUT, (Y(I), I=1,N)
C               TOUT = TOUT + 1.
C               IF (TOUT .LE. 10.) GO TO 20
C          100  FORMAT(...)
C               END (Sample)
C
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  SDRIV3,R1MACH,XERROR
C***END PROLOGUE  SDRIV2
      EXTERNAL F, JACOBN, FA, G
      REAL EPS, EWT, EWTCOM(1), G, HMAX, R1MACH, T, TOUT,
     8     WORK(*), Y(*)
      INTEGER IWORK(*)
      CHARACTER MSG*81
      PARAMETER(IMPL = 0, ML = 0, MU = 0, NDE = 0, MXSTEP = 1000)
C***FIRST EXECUTABLE STATEMENT  SDRIV2
      IF (MINT .LT. 1 .OR. MINT .GT. 3) THEN
        WRITE(MSG, '(''SDRIV21FE Illegal input.  Improper value for '',
     8  ''the integration method flag,'', I8)') MINT
        CALL XERROR(MSG, 81, 21, 2)
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
      IF (EWT .NE. 0.E0) THEN
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
      HMAX = SQRT(R1MACH(2))
      CALL SDRIV3 (N, T, Y, F, NSTATE, TOUT, NTASK, NROOT, EPS, EWTCOM,
     8             IERROR, MINT, MITER, IMPL, ML, MU, MXORD, HMAX, WORK,
     8             LENW, IWORK, LENIW, JACOBN, FA, NDE, MXSTEP, G)
      IF (MSTATE .GE. 0) THEN
        MSTATE = NSTATE
      ELSE
        MSTATE = - NSTATE
      END IF
      END
      SUBROUTINE SDRIV3 (N,T,Y,F,NSTATE,TOUT,NTASK,NROOT,EPS,EWT,IERROR,
     8   MINT,MITER,IMPL,ML,MU,MXORD,HMAX,WORK,LENW,IWORK,LENIW,JACOBN,
     8   FA,NDE,MXSTEP,G)
C***BEGIN PROLOGUE  SDRIV3
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850924   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***KEYWORDS  ODE,STIFF,ORDINARY DIFFERENTIAL EQUATIONS,
C             INITIAL VALUE PROBLEMS,GEAR'S METHOD,
C             SINGLE PRECISION
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***PURPOSE  The function of SDRIV3 is to solve N ordinary differential
C            equations of the form dY(I)/dT = F(Y(I),T), given the
C            initial conditions Y(I) = YI.  The program has options to
C            allow the solution of both stiff and non-stiff differential
C            equations.  Other important options are available.  SDRIV3
C            uses single precision arithmetic.
C***DESCRIPTION
C
C  I.  ABSTRACT  .......................................................
C
C    The primary function of SDRIV3 is to solve N ordinary differential
C    equations of the form dY(I)/dT = F(Y(I),T), given the initial
C    conditions Y(I) = YI.  The program has options to allow the
C    solution of both stiff and non-stiff differential equations.  In
C    addition, SDRIV3 may be used to solve:
C      1. The initial value problem, A*dY(I)/dT = F(Y(I),T), where A is
C         a non-singular matrix depending on Y and T.
C      2. The hybrid differential/algebraic initial value problem,
C         A*dY(I)/dT = F(Y(I),T), where A is a vector (whose values may
C         depend upon Y and T) some of whose components will be zero
C         corresponding to those equations which are algebraic rather
C         than differential.
C    SDRIV3 is to be called once for each output point of T.
C
C  II.  PARAMETERS  ....................................................
C
C    The user should use parameter names in the call sequence of SDRIV3
C    for those quantities whose value may be altered by SDRIV3.  The
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
C             JACOBN, FA, USERS, and G.
C
C    F      = A subroutine supplied by the user.  The name must be
C             declared EXTERNAL in the user's calling program.  This
C             subroutine is of the form:
C                   SUBROUTINE F (N, T, Y, YDOT)
C                   REAL Y(*), YDOT(*)
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
C             user's declaration in the program which calls SDRIV3.
C             Thus the dimensioning of Y in F, while required by FORTRAN
C             convention, does not actually allocate any storage.  When
C             this subroutine is called, the first N components of Y are
C             intermediate approximations to the solution components.
C             The user should not alter these values.  Here YDOT is a
C             vector of length N.  The user should only compute YDOT(I)
C             for I from 1 to N.
C
C    NSTATE = An integer describing the status of integration.  The
C             meaning of NSTATE is as follows:
C               1  (Input) Means the first call to the routine.  This
C                  value must be set by the user.  On all subsequent
C                  calls the value of NSTATE should be tested by the
C                  user, but must not be altered.  (As a convenience to
C                  the user who may wish to put out the initial
C                  conditions, SDRIV3 can be called with NSTATE=1, and
C                  TOUT=T.  In this case the program will return with
C                  NSTATE unchanged, i.e., NSTATE=1.)
C               2  (Output) Means a successful integration.  If a normal
C                  continuation is desired (i.e., a further integration
C                  in the same direction), simply advance TOUT and call
C                  again.  All other parameters are automatically set.
C               3  (Output)(Unsuccessful) Means the integrator has taken
C                  MXSTEP steps without reaching TOUT.  The user can
C                  continue the integration by simply calling SDRIV3
C                  again.
C               4  (Output)(Unsuccessful) Means too much accuracy has
C                  been requested.  EPS has been increased to a value
C                  the program estimates is appropriate.  The user can
C                  continue the integration by simply calling SDRIV3
C                  again.
C               5  (Output) A root was found at a point less than TOUT.
C                  The user can continue the integration toward TOUT by
C                  simply calling SDRIV3 again.
C
C    TOUT   = (Input) The point at which the solution is desired.  The
C             position of TOUT relative to T on the first call
C             determines the direction of integration.
C
C    NTASK  = (Input) An index specifying the manner of returning the
C             solution, according to the following:
C               NTASK = 1  Means SDRIV3 will integrate past TOUT and
C                          interpolate the solution.  This is the most
C                          efficient mode.
C               NTASK = 2  Means SDRIV3 will return the solution after
C                          each internal integration step, or at TOUT,
C                          whichever comes first.  In the latter case,
C                          the program integrates exactly to TOUT.
C               NTASK = 3  Means SDRIV3 will adjust its internal step to
C                          reach TOUT exactly (useful if a singularity
C                          exists beyond TOUT.)
C
C    NROOT  = (Input) The number of equations whose roots are desired.
C             If NROOT is zero, the root search is not active.  This
C             option is useful for obtaining output at points which are
C             not known in advance, but depend upon the solution, e.g.,
C             when some solution component takes on a specified value.
C             The root search is carried out using the user-written
C             function G (see description of G below.)  SDRIV3 attempts
C             to find the value of T at which one of the equations
C             changes sign.  SDRIV3 can find at most one root per
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
C             because the amount of work done by SDRIV3 increases as EPS
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
C                          by the user-written routine named USERS.
C                          This option allows all matrix algebra and
C                          storage decisions to be made by the user.
C                          The routine USERS is called by SDRIV3 when
C                          certain linear systems must be solved.  The
C                          user may choose any method to form, store and
C                          solve these systems in order to obtain the
C                          solution result that is returned to SDRIV3.
C                          In particular, this allows sparse matrix
C                          methods to be used.
C                          The call sequence for this routine is
C
C                           SUBROUTINE USERS (Y, YH, YWT, SAVE1, SAVE2,
C                          8              T, H, EL, IMPL, N, NDE, IFLAG)
C                           REAL Y(*), YH(*), YWT(*),
C                          8        SAVE1(*), SAVE2(*), T, H, EL
C
C                          The input variable IFLAG indicates what
C                          action is to be taken.  Subroutine USERS
C                          should perform the following operations,
C                          depending on the value of IFLAG and IMPL.
C
C                          IFLAG = 0
C                            IMPL = 0.  USERS is not called.
C                            IMPL = 1 or 2.  Solve the system
C                                A*X = SAVE2,
C                              returning the result in SAVE2.  The array
C                              SAVE1 can be used as a work array.
C
C                          IFLAG = 1
C                            IMPL = 0.  Compute, decompose and store the
C                              matrix (I - H*EL*J), where I is the
C                              identity matrix and J is the Jacobian
C                              matrix of the right hand side.  The array
C                              SAVE1 can be used as a work array.
C                            IMPL = 1 or 2. Compute, decompose and store
C                              the matrix (A - H*EL*J).  The array SAVE1
C                              can be used as a work array.
C
C                          IFLAG = 2
C                            IMPL = 0.   Solve the system
C                                (I - H*EL*J)*X = H*SAVE2 - YH - SAVE1,
C                              returning the result in SAVE2.
C                            IMPL = 1, or 2.  Solve the system
C                              (A - H*EL*J)*X = H*SAVE2 - A*(YH + SAVE1)
C                              returning the result in SAVE2.
C                            The array SAVE1 should not be altered.
C
C                          When using a value of MITER = 3, the
C                          subroutine FA is not required, even if IMPL
C                          is not 0.  For further information on using
C                          this option, see section IV-F below.
C
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
C             WORK is an array of LENW real words used
C             internally for temporary storage.  The user must allocate
C             space for this array in the calling program by a statement
C             such as
C                       REAL WORK(...)
C             The following table gives the required minimum value for
C             the length of WORK, depending on the value of IMPL and
C             MITER.  LENW should be set to the value used.  The
C             contents of WORK should not be disturbed between calls to
C             SDRIV3.
C
C      IMPL =   0                   1                   2
C              ---------------------------------------------------------
C MITER =  0   (MXORD+4)*N +       Not allowed         Not allowed
C              2*NROOT + 204
C
C         1,2  N*N+(MXORD+4)*N     2*N*N+(MXORD+4)*N   N*N+(MXORD+5)*N
C              + 2*NROOT + 204     + 2*NROOT + 204     + 2*NROOT + 204
C
C          3   (MXORD+4)*N +       (MXORD+4)*N +       (MXORD+4)*N +
C              2*NROOT + 204       2*NROOT + 204       2*NROOT + 204
C
C         4,5  (2*ML+MU)*N +       (4*ML+2*MU)*N +     (2*ML+MU)*N +
C              (MXORD+5)*N +       (MXORD+6)*N +       (MXORD+6)*N +
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
C             of IWORK should not be disturbed between calls to SDRIV3.
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
C             approximated numerically by SDRIV3, and we only ask the
C             user to tell SDRIV3 the value of ML and MU if the system
C             is banded.  For the cases MITER = 1 and 4 the user must
C             derive these partials algebraically and encode them in
C             subroutine JACOBN.  By computing these derivatives the
C             user can often save 20-30 per cent of the computing time.
C             Usually, however, the accuracy is not much affected and
C             most users will probably forego this option.  The optional
C             user-written subroutine JACOBN has the form:
C                   SUBROUTINE JACOBN (N, T, Y, DFDY, MATDIM, ML, MU)
C                   REAL Y(*), DFDY(MATDIM,*)
C                     .
C                     .
C                     Calculate values of DFDY
C                     .
C                     .
C                   END (Sample)
C             Here Y is a vector of length at least N.  The actual
C             length of Y is determined by the user's declaration in the
C             program which calls SDRIV3.  Thus the dimensioning of Y in
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
C             K=I-J+MU+1.
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
C                   REAL Y(*), A(MATDIM,*)
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
C                   REAL Y(*), A(*)
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
C             program which calls SDRIV3.  Thus the dimensioning of Y in
C             FA, while required by FORTRAN convention, does not
C             actually allocate any storage.  When this subroutine is
C             called, the first N components of Y are intermediate
C             approximations to the solution components.  The user
C             should not alter these values.  FA is always called
C             immediately after calling F, with the same values of T
C             and Y.
C
C    NDE    = (Input) The number of differential equations.  This is
C             required only for IMPL = 2, with NDE .LT. N.
C
C    MXSTEP = (Input) The maximum number of internal steps allowed on
C             one call to SDRIV3.
C
C    G      = A real FORTRAN function supplied by the user
C             if NROOT is not 0.  In this case, the name must be
C             declared EXTERNAL in the user's calling program.  G is
C             repeatedly called with different values of IROOT to obtain
C             the value of each of the NROOT equations for which a root
C             is desired.  G is of the form:
C                   REAL FUNCTION G (N, T, Y, IROOT)
C                   REAL Y(*)
C                   GO TO (10, ...), IROOT
C              10   G = ...
C                     .
C                     .
C                   END (Sample)
C             Here, Y is a vector of length at least N, whose first N
C             components are the solution components at the point T.
C             The user should not alter these values.  The actual length
C             of Y is determined by the user's declaration in the
C             program which calls SDRIV3.  Thus the dimensioning of Y in
C             G, while required by FORTRAN convention, does not actually
C             allocate any storage.
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
C       come from SDRIV3:
C
C        No.  Type         Message
C        ---  ----         -------
C         1   Fatal        From SDRIV2: The integration method flag has
C                          an illegal value.
C         2   Warning      The output point is inconsistent with the
C                          value of NTASK and T.
C         3   Warning      Number of steps to reach TOUT exceeds MXSTEP.
C         4   Recoverable  Requested accuracy is too stringent.
C         5   Warning      Step size is below the roundoff level.
C         6   Fatal        EPS is less than zero.
C         7   Fatal        N is not positive.
C         8   Fatal        Insufficient work space provided.
C         9   Fatal        Improper value for MINT, MITER and/or IMPL.
C        10   Fatal        The IWORK array is too small.
C        11   Fatal        The step size has gone to zero.
C        12   Fatal        Excessive amount of work.
C        13   Fatal        For IMPL=1 or 2, the matrix A is singular.
C        14   Fatal        MXORD is not positive.
C        15   Fatal        From SDRIV1: N is greater than 200.
C        16   Fatal        From SDRIV1: The WORK array is too small.
C
C    B. The first three elements of WORK and the first five elements of
C       IWORK will contain the following statistical data:
C         AVGH     The average step size used.
C         HUSED    The step size last used (successfully).
C         AVGORD   The average order used.
C         IMXERR   The index of the element of the solution vector that
C                  contributed most to the last error test.
C         NQUSED   The order last used (successfully).
C         NSTEP    The number of steps taken.
C         NFE      The number of evaluations of the right hand side.
C         NJE      The number of evaluations of the Jacobian matrix.
C
C  IV.  REMARKS  .......................................................
C
C    A. Other routines used:
C         SDNTP, SDZRO, SDSTP, SDNTL, SDPST, SDCOR, SDCST,
C         SDPSC, and SDSCL;
C         SGEFA, SGESL, SGBFA, SGBSL, and SNRM2 (from LINPACK)
C         R1MACH (from the Bell Laboratories Machine Constants Package)
C         XERROR (from the SLATEC Common Math Library)
C       The last seven routines above, not having been written by the
C       present authors, are not explicitly part of this package.
C
C    B. On any return from SDRIV3 all information necessary to continue
C       the calculation is contained in the call sequence parameters,
C       including the work arrays.  Thus it is possible to suspend one
C       problem, integrate another, and then return to the first.
C
C    C. There are user-written routines which are only required by
C       SDRIV3 when certain parameters are set.  Thus a message warning
C       of unsatisfied externals may be issued during the load or link
C       phase.  This message should never refer to F.  This message can
C       be ignored if: it refers to JACOBN and MITER is not 1 or 4, or
C       it refers to FA and IMPL is 0 or MITER is 3, or it refers to
C       USERS and MITER is not 3, or it refers to G and NROOT is 0.
C
C    D. If this package is to be used in an overlay situation, the user
C       must declare in the primary overlay the variables in the call
C       sequence to SDRIV3.
C
C    E. Changing parameters during an integration.
C       The value of NROOT, EPS, EWT, IERROR, MINT, MITER, or HMAX may
C       be altered by the user between calls to SDRIV3.  For example, if
C       too much accuracy has been requested (the program returns with
C       NSTATE = 4 and an increased value of EPS) the user may wish to
C       increase EPS further.  In general, prudence is necessary when
C       making changes in parameters since such changes are not
C       implemented until the next integration step, which is not
C       necessarily the next call to SDRIV3.  This can happen if the
C       program has already integrated to a point which is beyond the
C       new point TOUT.
C
C    F. As the price for complete control of matrix algebra, the SDRIV3
C       USERS option puts all responsibility for Jacobian matrix
C       evaluation on the user.  It is often useful to approximate
C       numerically all or part of the Jacobian matrix.  However this
C       must be done carefully.  The FORTRAN sequence below illustrates
C       the method we recommend.  It can be inserted directly into
C       subroutine USERS to approximate Jacobian elements in rows I1
C       to I2 and columns J1 to J2.
C              REAL DFDY(N,N), EPSJ, H, R, R1MACH,
C             8     SAVE1(N), SAVE2(N), T, UROUND, Y(N), YJ, YWT(N)
C              UROUND = R1MACH(4)
C              EPSJ = UROUND**(1.E0/3.E0)
C              DO 30 J = J1,J2
C                R = EPSJ*MAX(ABS(YWT(J)), ABS(Y(J)))
C                IF (R .EQ. 0.E0) R = EPSJ
C                YJ = Y(J)
C                Y(J) = Y(J) + R
C                CALL F (N, T, Y, SAVE1)
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
C***REFERENCES  GEAR, C. W., "NUMERICAL INITIAL VALUE PROBLEMS IN
C                 ORDINARY DIFFERENTIAL EQUATIONS", PRENTICE-HALL, 1971.
C***ROUTINES CALLED  SDSTP,SDNTP,SDZRO,SGEFA,SGESL,SGBFA,SGBSL,SNRM2,
C                    R1MACH,XERROR
C***END PROLOGUE  SDRIV3
      EXTERNAL F, JACOBN, FA, G
      REAL AE, BIG, EPS, EWT(*), G, GLAST, H, HMAX, HSIGN,
     8     NROUND, RE, R1MACH, SIZE, SNRM2, SUM, T, TLAST, TOUT, TROOT,
     8     UROUND, WORK(*), Y(*)
      INTEGER IWORK(*)
      LOGICAL CONVRG
      CHARACTER MSG*205
      PARAMETER(NROUND = 20.E0)
      PARAMETER(IAVGH = 1, IHUSED = 2, IAVGRD = 3,
     8          IEL = 4, IH = 160, IHMAX = 161, IHOLD = 162,
     8          IHSIGN = 163, IRC = 164, IRMAX = 165, IT = 166,
     8          ITOUT = 167, ITQ = 168, ITREND = 204, IYH = 205,
     8          INDMXR = 1, INQUSD = 2, INSTEP = 3, INFE = 4, INJE = 5,
     8          INROOT = 6, ICNVRG = 7, IJROOT = 8, IJTASK = 9,
     8          IMNTLD = 10, IMTRLD = 11, INQ = 12, INRTLD = 13,
     8          INDTRT = 14, INWAIT = 15, IMNT = 16, IMTRSV = 17,
     8          IMTR = 18, IMXRDS = 19, IMXORD = 20)
      PARAMETER(INDPRT = 21, INDPVT = 22)
C***FIRST EXECUTABLE STATEMENT  SDRIV3
      UROUND = R1MACH (4)
      IF (NROOT .NE. 0) THEN
        AE = R1MACH(1)
        RE = UROUND
      END IF
      IF (EPS .LT. 0.E0) THEN
        WRITE(MSG, '(''SDRIV36FE Illegal input.  EPS,'', E16.8,
     8  '', is negative.'')') EPS
        CALL XERROR(MSG, 60, 6, 2)
        RETURN
      END IF
      IF (N .LE. 0) THEN
        WRITE(MSG, '(''SDRIV37FE Illegal input.  Number of equations,'',
     8  I8, '', is not positive.'')') N
        CALL XERROR(MSG, 72, 7, 2)
        RETURN
      END IF
      IF (MXORD .LE. 0) THEN
        WRITE(MSG, '(''SDRIV314FE Illegal input.  Maximum order,'', I8,
     8  '', is not positive.'')') MXORD
        CALL XERROR(MSG, 67, 14, 2)
        RETURN
      END IF
      IF ((MINT .LT. 1 .OR. MINT .GT. 3) .OR. (MINT .EQ. 3 .AND.
     8  (MITER .EQ. 0 .OR. MITER .EQ. 3 .OR. IMPL .NE. 0))
     8  .OR. (MITER .LT. 0 .OR. MITER .GT. 5) .OR.
     8  (IMPL .NE. 0 .AND. IMPL .NE. 1 .AND. IMPL .NE. 2) .OR.
     8  ((IMPL .EQ. 1 .OR. IMPL .EQ. 2) .AND. MITER .EQ. 0) .OR.
     8  (IMPL .EQ. 2 .AND. MINT .EQ. 1)) THEN
        WRITE(MSG, '(''SDRIV39FE Illegal input.  Improper value for '',
     8  ''MINT, MITER and/or IMPL.'')')
        CALL XERROR(MSG, 69, 9, 2)
        RETURN
      END IF
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) THEN
        LIWCHK = INDPVT - 1
      ELSE IF (MITER .EQ. 1 .OR. MITER .EQ. 2 .OR. MITER .EQ. 4 .OR.
     8  MITER .EQ. 5) THEN
        LIWCHK = INDPVT + N - 1
      END IF
      IF (LENIW .LT. LIWCHK) THEN
        WRITE(MSG, '(''SDRIV310FE Illegal input.  Insufficient '',
     8  ''storage allocated for the IWORK array.  Based on the '')')
        WRITE(MSG(94:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LIWCHK
        CALL XERROR(MSG, 164, 10, 2)
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
      IA = ITROOT + NROOT
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
        WRITE(MSG, '(''SDRIV38FE Illegal input.  Insufficient '',
     8  ''storage allocated for the WORK array.  Based on the '')')
        WRITE(MSG(92:), '(''value of the input parameters involved, '',
     8  ''the required storage is'', I8)') LENCHK
        CALL XERROR(MSG, 162, 8, 2)
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
        IF (T .EQ. TOUT) RETURN
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
        H = (TOUT - T)*(1.E0 - 4.E0*UROUND)
        H = SIGN(MIN(ABS(H), HMAX), H)
        WORK(IH) = H
        HSIGN = SIGN(1.E0, H)
        WORK(IHSIGN) = HSIGN
        IWORK(IJTASK) = 0
        WORK(IAVGH) = 0.E0
        WORK(IAVGRD) = 0.E0
        IWORK(INQUSD) = 0
        IWORK(INSTEP) = 0
        IWORK(INFE) = 0
        IWORK(INJE) = 0
        WORK(IT) = T
        IWORK(ICNVRG) = 0
        IWORK(INDPRT) = 0
C                                                 Set initial conditions
        DO 30 I = 1,N
          JYH = I + IYH - 1
 30       WORK(JYH) = Y(I)
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
      IF(NROOT .NE. 0) THEN
        JROOT = IWORK(IJROOT)
        IF (JROOT .GT. 0) THEN
C                                      TOUT has just been reported.
C                                      If TROOT .LE. TOUT, report TROOT.
          IF (NSTATE .NE. 5) THEN
            IF (TOUT*HSIGN .GE. WORK(ITOUT)*HSIGN) THEN
              TROOT = WORK(ITOUT)
              CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
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
                  CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT,WORK(IYH),Y)
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
                CALL SDNTP(H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH), Y)
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
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
      ELSE IF (NTASK .EQ. 2) THEN
C                                                      Check if TOUT has
C                                                      been reset .LT. T
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''SDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
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
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.E0) GO TO 670
          IWORK(IJTASK) = -1
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        NSTATE = 2
        IF (T*HSIGN .GT. TOUT*HSIGN) THEN
          WRITE(MSG, '(''SDRIV32WRN With NTASK='', I1, '' on input, '',
     8    ''T,'', E16.8, '', was beyond TOUT,'', E16.8, ''.  Solution'',
     8    '' obtained by interpolation.'')') NTASK, T, TOUT
          CALL XERROR(MSG, 124, 2, 0)
          CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
          T = TOUT
          GO TO 580
        END IF
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
          GO TO 560
        END IF
        IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
          H = TOUT - T
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
          WORK(IH) = H
          IF (H .EQ. 0.E0) GO TO 670
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
 200      WORK(JGNOW) = G (N, T, Y, I)
      END IF
      IF (IERROR .EQ. 1) THEN
        DO 230 I = 1,N
          JYWT = I + IYWT - 1
 230      WORK(JYWT) = 1.E0
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
          IF (Y(I) .EQ. 0.E0) GO TO 290
          JYWT = I + IYWT - 1
 280      WORK(JYWT) = ABS(Y(I))
        GO TO 410
 290    IF (IWORK(IJTASK) .EQ. 0) THEN
          CALL F (N, T, Y, WORK(ISAVE2))
          IWORK(INFE) = IWORK(INFE) + 1
          IF (MITER .EQ. 3 .AND. IMPL .NE. 0) THEN
            IFLAG = 0
            CALL USERS(Y,WORK(IYH),WORK(IYWT),WORK(ISAVE1),WORK(ISAVE2),
     8                 T, H, WORK(IEL), IMPL, N, NDECOM, IFLAG)
          ELSE IF (IMPL .EQ. 1) THEN
            IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
              CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
              CALL SGEFA (WORK(IA), MATDIM, N, IWORK(INDPVT), INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL SGESL(WORK(IA),MATDIM,N,IWORK(INDPVT),WORK(ISAVE2),0)
            ELSE IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
              JAML = IA + ML
              CALL FA (N, T, Y, WORK(JAML), MATDIM, ML, MU, NDECOM)
              CALL SGBFA (WORK(IA),MATDIM,N,ML,MU,IWORK(INDPVT),INFO)
              IF (INFO .NE. 0) GO TO 690
              CALL SGBSL (WORK(IA), MATDIM, N, ML, MU, IWORK(INDPVT),
     8                    WORK(ISAVE2), 0)
            END IF
          ELSE IF (IMPL .EQ. 2) THEN
            CALL FA (N, T, Y, WORK(IA), MATDIM, ML, MU, NDECOM)
            DO 340 I = 1,NDECOM
              JA = I + IA - 1
              JSAVE2 = I + ISAVE2 - 1
              IF(WORK(JA) .EQ. 0.E0) GO TO 690
 340          WORK(JSAVE2) = WORK(JSAVE2)/WORK(JA)
          END IF
        END IF
        DO 360 J = I,N
          JYWT = J + IYWT - 1
          IF (Y(J) .NE. 0.E0) THEN
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
          IF (WORK(JYWT) .EQ. 0.E0) WORK(JYWT) = UROUND
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
      SUM = SNRM2(N, WORK(ISAVE2), 1)/SQRT(REAL(N))
      IF (EPS .LT. SUM*UROUND) THEN
        EPS = SUM*UROUND*(1.E0 + 10.E0*UROUND)
        WRITE(MSG, '(''SDRIV34REC At T,'', E16.8, '', the requested '',
     8  ''accuracy, EPS, was not obtainable with the machine '',
     8  ''precision.  EPS has been increased to'')') T
        WRITE(MSG(137:), '(E16.8)') EPS
        CALL XERROR(MSG, 152, 4, 1)
        NSTATE = 4
        GO TO 560
      END IF
      IF (ABS(H) .GE. UROUND*ABS(T)) THEN
        IWORK(INDPRT) = 0
      ELSE IF (IWORK(INDPRT) .EQ. 0) THEN
        WRITE(MSG, '(''SDRIV35WRN At T,'', E16.8, '', the step size,'',
     8  E16.8, '', is smaller than the roundoff level of T.  '')') T, H
        WRITE(MSG(109:), '(''This may occur if there is an abrupt '',
     8  ''change in the right hand side of the differential '',
     8  ''equations.'')')
        CALL XERROR(MSG, 205, 5, 0)
        IWORK(INDPRT) = 1
      END IF
      IF (NTASK.NE.2) THEN
        IF ((IWORK(INSTEP)-NSTEPL) .GT. MXSTEP) THEN
          WRITE(MSG, '(''SDRIV33WRN At T,'', E16.8, '', '', I8,
     8    '' steps have been taken without reaching TOUT,'', E16.8)')
     8    T, MXSTEP, TOUT
          CALL XERROR(MSG, 103, 3, 0)
          NSTATE = 3
          GO TO 560
        END IF
      END IF
C
C     CALL SDSTP (EPS, F, FA, HMAX, IMPL, JACOBN, MATDIM, MAXORD,
C    8            MINT, MITER, ML, MU, N, NDE, YWT, UROUND,
C    8            AVGH, AVGORD, H, HUSED, JTASK, MNTOLD, MTROLD,
C    8            NFE, NJE, NQUSED, NSTEP, T, Y, YH,  A, CONVRG,
C    8            DFDY, EL, HOLD, IPVT, JSTATE, NQ, NWAIT, RC,
C    8            RMAX, SAVE1, SAVE2, TQ, TREND, ISWFLG, MTRSV, MXRDSV)
C
      CALL SDSTP (EPS, F, FA, WORK(IHMAX), IMPL, JACOBN, MATDIM,
     8            IWORK(IMXORD), IWORK(IMNT), IWORK(IMTR), ML, MU, N,
     8            NDECOM, WORK(IYWT), UROUND, WORK(IAVGH), WORK(IAVGRD),
     8            WORK(IH), WORK(IHUSED), IWORK(IJTASK), IWORK(IMNTLD),
     8            IWORK(IMTRLD), IWORK(INFE), IWORK(INJE),
     8            IWORK(INQUSD), IWORK(INSTEP), WORK(IT), Y, WORK(IYH),
     8            WORK(IA), CONVRG, WORK(IDFDY), WORK(IEL), WORK(IHOLD),
     8            IWORK(INDPVT), JSTATE, IWORK(INQ), IWORK(INWAIT),
     8            WORK(IRC), WORK(IRMAX), WORK(ISAVE1), WORK(ISAVE2),
     8            WORK(ITQ), WORK(ITREND), MINT, IWORK(IMTRSV),
     8            IWORK(IMXRDS))
      T = WORK(IT)
      H = WORK(IH)
      GO TO (470, 670, 680, 690), JSTATE
 470  IWORK(IJTASK) = 1
C                                 Determine if a root has been overtaken
      IF (ABS(H) .GE. UROUND*ABS(T) .AND. NROOT .NE. 0) THEN
        IROOT = 0
        DO 500 I = 1,NROOT
          JTROOT = ITROOT + I - 1
          JGNOW = IGNOW + I - 1
          GLAST = WORK(JGNOW)
          WORK(JGNOW) = G (N, T, Y, I)
          IF (GLAST*WORK(JGNOW) .GT. 0.E0) THEN
            WORK(JTROOT) = T + H
          ELSE
            IF (WORK(JGNOW) .EQ. 0.E0) THEN
              WORK(JTROOT) = T
              IROOT = I
            ELSE
              IF (GLAST .EQ. 0.E0) THEN
                WORK(JTROOT) = T + H
              ELSE
                TLAST = T - WORK(IHUSED)
                IROOT = I
                TROOT = T
                CALL SDZRO (AE, G, H, N, IWORK(INQ), IROOT, RE, T,
     8                      WORK(IYH), UROUND,  TROOT, TLAST,
     8                      WORK(JGNOW), GLAST,  WORK(ISAVE1))
                WORK(JTROOT) = TROOT
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
            CALL SDNTP (H, 0, N, IWORK(INQ), T, TROOT, WORK(IYH),  Y)
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
        CALL SDNTP (H, 0, N, IWORK(INQ), T, TOUT, WORK(IYH),  Y)
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
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.E0) GO TO 670
            IWORK(IJTASK) = -1
          END IF
        END IF
      ELSE IF (NTASK .EQ. 3) THEN
        IF (ABS(TOUT - T).LE.NROUND*UROUND*MAX(ABS(T), ABS(TOUT))) THEN
          T = TOUT
        ELSE
          IF ((T + H)*HSIGN .GT. TOUT*HSIGN) THEN
            H = TOUT - T
            IF ((T + H)*HSIGN.GT.TOUT*HSIGN) H = H*(1.E0 - 4.E0*UROUND)
            WORK(IH) = H
            IF (H .EQ. 0.E0) GO TO 670
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
      BIG = 0.E0
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
C                                        Fatal errors are processed here
C
 670  WRITE(MSG, '(''SDRIV311FE At T,'', E16.8, '', the attempted '',
     8  ''step size has gone to zero.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')') T
      CALL XERROR(MSG, 129, 11, 2)
      RETURN
C
 680  WRITE(MSG, '(''SDRIV312FE At T,'', E16.8, '', the step size has'',
     8  '' been reduced about 50 times without advancing the '')') T
      WRITE(MSG(103:), '(''solution.  Often this occurs if the '',
     8  ''problem setup is incorrect.'')')
      CALL XERROR(MSG, 165, 12, 2)
      RETURN
C
 690  WRITE(MSG, '(''SDRIV313FE At T,'', E16.8, '', while solving'',
     8  '' A*YDOT = F, A is singular.'')') T
      CALL XERROR(MSG, 74, 13, 2)
      RETURN
      END
      SUBROUTINE SDSCL (HMAX,N,NQ,RMAX,H,RC,RH,YH)
C***BEGIN PROLOGUE  SDSCL
C***REFER TO  SDRIV3
C   This subroutine rescales the YH array whenever the step size
C   is changed.
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  850319   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDSCL
      REAL H, HMAX, RC, RH, RMAX, R1, YH(N,*)
C***FIRST EXECUTABLE STATEMENT  SDSCL
      IF (H .LT. 1.E0) THEN
        RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
      ELSE
        RH = MIN(RH, RMAX, HMAX/ABS(H))
      END IF
      R1 = 1.E0
      DO 10 J = 1,NQ
        R1 = R1*RH
        DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
      H = H*RH
      RC = RC*RH
      END
      SUBROUTINE SDSTP (EPS,F,FA,HMAX,IMPL,JACOBN,MATDIM,MAXORD,MINT,
     8   MITER,ML,MU,N,NDE,YWT,UROUND,AVGH,AVGORD,H,HUSED,JTASK,MNTOLD,
     8   MTROLD,NFE,NJE,NQUSED,NSTEP,T,Y,YH,A,CONVRG,DFDY,EL,HOLD,IPVT,
     8   JSTATE,NQ,NWAIT,RC,RMAX,SAVE1,SAVE2,TQ,TREND,ISWFLG,MTRSV,
     8   MXRDSV)
C***BEGIN PROLOGUE  SDSTP
C***REFER TO  SDRIV3
C  SDSTP performs one step of the integration of an initial value
C  problem for a system of ordinary differential equations.
C  Communication with SDSTP is done with the following variables:
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
C              the value used in SDSTP.
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
C***ROUTINES CALLED  SDNTL,SDPST,SDCOR,SDPSC,SDSCL,SNRM2
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  860513   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDSTP
      EXTERNAL F, JACOBN, FA
      REAL A(MATDIM,*), AVGH, AVGORD, BIAS1, BIAS2, BIAS3,
     8     BND, CTEST, D, DENOM, DFDY(MATDIM,*), D1, EL(13,12), EPS,
     8     ERDN, ERUP, ETEST, H, HMAX, HN, HOLD, HS, HUSED, NUMER, RC,
     8     RCTEST, RH, RH1, RH2, RH3, RMAX, RMFAIL, RMNORM, SAVE1(*),
     8     SAVE2(*), SNRM2, T, TOLD, TQ(3,12), TREND, TRSHLD, UROUND,
     8     Y(*), YH(N,*), YWT(*), Y0NRM
      INTEGER IPVT(*)
      LOGICAL CONVRG, EVALFA, EVALJC, IER, SWITCH
      PARAMETER(BIAS1 = 1.3E0, BIAS2 = 1.2E0, BIAS3 = 1.4E0, MXFAIL = 3,
     8          MXITER = 3, MXTRY = 50, RCTEST = .3E0, RMFAIL = 2.E0,
     8          RMNORM = 10.E0, TRSHLD = 1.E0)
C***FIRST EXECUTABLE STATEMENT  SDSTP
      BND = 0.E0
      SWITCH = .FALSE.
      NTRY = 0
      TOLD = T
      NFAIL = 0
      IF (JTASK .LE. 0) THEN
        CALL SDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
      END IF
      IER = .FALSE.
 100  NTRY = NTRY + 1
      IF (NTRY .GT. MXTRY) GO TO 410
      T = T + H
      CALL SDPSC (1, N, NQ,  YH)
      EVALJC = ((ABS(RC - 1.E0) .GT. RCTEST) .AND. (MITER .NE. 0))
      EVALFA = .NOT. EVALJC
C
 110  ITER = 0
      DO 115 I = 1,N
 115    Y(I) = YH(I,1)
      CALL F (N, T, Y, SAVE2)
      NFE = NFE + 1
      IF (EVALJC .OR. IER) THEN
        CALL SDPST (EL, F, FA, H, IMPL, JACOBN, MATDIM, MITER, ML,
     8              MU, N, NDE, NQ, SAVE2, T, Y, YH, YWT, UROUND,
     8              NFE, NJE,  A, DFDY, IER, IPVT, SAVE1, ISWFLG, BND)
        IF (IER) GO TO 160
        CONVRG = .FALSE.
        RC = 1.E0
      END IF
      DO 125 I = 1,N
 125    SAVE1(I) = 0.E0
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
 130  CALL SDCOR (DFDY, EL, FA, H, IMPL, IPVT, MATDIM, MITER, ML,
     8            MU, N, NDE, NQ, T, Y, YH, YWT,  EVALFA, SAVE1,
     8            SAVE2,  A, D)
      IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
        IF (ITER .EQ. 0) THEN
          NUMER = SNRM2(N, SAVE1, 1)
          DO 132 I = 1,N
 132        DFDY(1,I) = SAVE1(I)
          Y0NRM = SNRM2(N, YH, 1)
        ELSE
          DENOM = NUMER
          DO 134 I = 1,N
 134        DFDY(1,I) = SAVE1(I) - DFDY(1,I)
          NUMER = SNRM2(N, DFDY, MATDIM)
          IF (EL(1,NQ)*NUMER .LE. 100.E0*UROUND*Y0NRM) THEN
            IF (RMAX .EQ. RMFAIL) THEN
              SWITCH = .TRUE.
              GO TO 170
            END IF
          END IF
          DO 136 I = 1,N
 136        DFDY(1,I) = SAVE1(I)
          IF (DENOM .NE. 0.E0)
     8    BND = MAX(BND, NUMER/(DENOM*ABS(H)*EL(1,NQ)))
        END IF
      END IF
      IF (ITER .GT. 0) TREND = MAX(.9E0*TREND, D/D1)
      D1 = D
      CTEST = MIN(2.E0*TREND, 1.E0)*D
      IF (CTEST .LE. EPS) GO TO 170
      ITER = ITER + 1
      IF (ITER .LT. MXITER) THEN
        DO 140 I = 1,N
 140      Y(I) = YH(I,1) + EL(1,NQ)*SAVE1(I)
        CALL F (N, T, Y, SAVE2)
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
      CALL SDPSC (-1, N, NQ,  YH)
      NWAIT = NQ + 2
      IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
      IF (ITER .EQ. 0) THEN
        RH = .3E0
      ELSE
        RH = .9E0*(EPS/CTEST)**(.2E0)
      END IF
      IF (RH*H .EQ. 0.E0) GO TO 400
      CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
      GO TO 100
C                          The corrector has converged.  CONVRG is set
C                          to .TRUE. if partial derivatives were used,
C                          to indicate that they may need updating on
C                          subsequent steps.  The error test is made.
 170  CONVRG = (MITER .NE. 0)
      DO 180 I = 1,NDE
 180    SAVE2(I) = SAVE1(I)/YWT(I)
      ETEST = SNRM2(NDE, SAVE2, 1)/(TQ(2,NQ)*SQRT(REAL(NDE)))
C
C                           The error test failed.  NFAIL keeps track of
C                           multiple failures.  Restore T and the YH
C                           array to their previous values, and prepare
C                           to try the step again.  Compute the optimum
C                           step size for this or one lower order.
      IF (ETEST .GT. EPS) THEN
        T = TOLD
        CALL SDPSC (-1, N, NQ,  YH)
        NFAIL = NFAIL + 1
        IF (NFAIL .LT. MXFAIL) THEN
          IF (JTASK .NE. 0 .AND. JTASK .NE. 2) RMAX = RMFAIL
          RH2 = 1.E0/(BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
          IF (NQ .GT. 1) THEN
            DO 190 I = 1,NDE
 190          SAVE2(I) = YH(I,NQ+1)/YWT(I)
            ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
            RH1 = 1.E0/MAX(1.E0, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
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
          IF (RH*H .EQ. 0.E0) GO TO 400
          CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
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
        CALL SDNTL (EPS, F, FA, HMAX, HOLD, IMPL, JTASK, MATDIM,
     8              MAXORD, MINT, MITER, ML, MU, N, NDE, SAVE1, T,
     8              Y, YWT,  H, MNTOLD, MTROLD, NFE, RC, YH,
     8              A, CONVRG, EL, IER, IPVT, NQ, NWAIT, RH, RMAX,
     8              SAVE2, TQ, TREND, ISWFLG)
        IF (H .EQ. 0.E0) GO TO 400
        IF (IER) GO TO 420
        GO TO 100
      END IF
C                          After a successful step, update the YH array.
      NSTEP = NSTEP + 1
      HUSED = H
      NQUSED = NQ
      AVGH = (REAL(NSTEP-1)*AVGH + H)/REAL(NSTEP)
      AVGORD = (REAL(NSTEP-1)*AVGORD + REAL(NQ))/REAL(NSTEP)
      DO 230 J = 1,NQ+1
        DO 230 I = 1,N
 230      YH(I,J) = YH(I,J) + EL(J,NQ)*SAVE1(I)
      DO 235 I = 1,N
 235    Y(I) = YH(I,1)
C                                          If ISWFLG is 3, consider
C                                          changing integration methods.
C
      IF (ISWFLG .EQ. 3) THEN
        IF (BND .NE. 0.E0) THEN
          IF (MINT .EQ. 1 .AND. NQ .LE. 5) THEN
            HN = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            HS = ABS(H)/MAX(UROUND,
     8      (ETEST/(EPS*EL(NQ+1,1)))**(1.E0/REAL(NQ+1)))
            IF (HS .GT. 1.2E0*HN) THEN
              MINT = 2
              MNTOLD = MINT
              MITER = MTRSV
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 5)
              RC = 0.E0
              RMAX = RMNORM
              TREND = 1.E0
              CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
              NWAIT = NQ + 2
            END IF
          ELSE IF (MINT .EQ. 2) THEN
            HS = ABS(H)/MAX(UROUND, (ETEST/EPS)**(1.E0/REAL(NQ+1)))
            HN = ABS(H)/MAX(UROUND,
     8      (ETEST*EL(NQ+1,1)/EPS)**(1.E0/REAL(NQ+1)))
            HN = MIN(HN, 1.E0/(2.E0*EL(1,NQ)*BND))
            IF (HN .GE. HS) THEN
              MINT = 1
              MNTOLD = MINT
              MITER = 0
              MTROLD = MITER
              MAXORD = MIN(MXRDSV, 12)
              RMAX = RMNORM
              TREND = 1.E0
              CONVRG = .FALSE.
              CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
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
        RC = 0.E0
        RMAX = RMNORM
        TREND = 1.E0
        CALL SDCST (MAXORD, MINT, ISWFLG, EL, TQ)
        NWAIT = NQ + 2
      END IF
C                           Consider changing H if NWAIT = 1.  Otherwise
C                           decrease NWAIT by 1.  If NWAIT is then 1 and
C                           NQ.LT.MAXORD, then SAVE1 is saved for use in
C                           a possible order increase on the next step.
C
      IF (JTASK .EQ. 0 .OR. JTASK .EQ. 2) THEN
        RH = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (RH.GT.TRSHLD) CALL SDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
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
          RH1 = 0.E0
        ELSE
          DO 270 I = 1,NDE
 270        SAVE2(I) = YH(I,NQ+1)/YWT(I)
          ERDN = SNRM2(NDE, SAVE2, 1)/(TQ(1,NQ)*SQRT(REAL(NDE)))
          RH1 = 1.E0/MAX(UROUND, BIAS1*(ERDN/EPS)**(1.E0/REAL(NQ)))
        END IF
        RH2 = 1.E0/MAX(UROUND, BIAS2*(ETEST/EPS)**(1.E0/REAL(NQ+1)))
        IF (NQ .EQ. MAXORD) THEN
          RH3 = 0.E0
        ELSE
          DO 290 I = 1,NDE
 290        SAVE2(I) = (SAVE1(I) - YH(I,MAXORD+1))/YWT(I)
          ERUP = SNRM2(NDE, SAVE2, 1)/(TQ(3,NQ)*SQRT(REAL(NDE)))
          RH3 = 1.E0/MAX(UROUND, BIAS3*(ERUP/EPS)**(1.E0/REAL(NQ+2)))
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
 360        YH(I,NQ+2) = SAVE1(I)*EL(NQ+1,NQ)/REAL(NQ+1)
          NQ = NQ + 1
          RC = RC*EL(1,NQ)/EL(1,NQ-1)
        END IF
        IF (ISWFLG .EQ. 3 .AND. MINT .EQ. 1) THEN
          IF (BND.NE.0.E0) RH = MIN(RH, 1.E0/(2.E0*EL(1,NQ)*BND*ABS(H)))
        END IF
        CALL SDSCL (HMAX, N, NQ, RMAX,  H, RC, RH, YH)
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
      END
      SUBROUTINE SDZRO (AE,F,H,N,NQ,IROOT,RE,T,YH,UROUND,B,C,FB,FC,Y)
C***BEGIN PROLOGUE  SDZRO
C***REFER TO  SDRIV3
C     This is a special purpose version of ZEROIN, modified for use with
C     the SDRIV1 package.
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
C                real result.  This name must be in an
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
C***ROUTINES CALLED  SDNTP
C***DATE WRITTEN   790601   (YYMMDD)
C***REVISION DATE  841119   (YYMMDD)
C***CATEGORY NO.  I1A2,I1A1B
C***AUTHOR  KAHANER, D. K., NATIONAL BUREAU OF STANDARDS,
C           SUTHERLAND, C. D., LOS ALAMOS NATIONAL LABORATORY
C***END PROLOGUE  SDZRO
      EXTERNAL F
      REAL A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC,
     8     H, P, Q, RE, RW, T, TOL, UROUND, Y(*), YH(N,*)
C***FIRST EXECUTABLE STATEMENT  SDZRO
      ER = 4.E0*UROUND
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
      CMB = 0.5E0*(C - B)
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
      IF (P .LT. 0.E0) THEN
        P = -P
        Q = -Q
      END IF
C                          Update A and check for satisfactory reduction
C                          in the size of our bounding interval.
      A = B
      FA = FB
      IC = IC + 1
      IF (IC .GE. 4) THEN
        IF (8.E0*ACMB .GE. ACBS) THEN
C                                                                 Bisect
          B = 0.5E0*(C + B)
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
        B = 0.5E0*(C + B)
      END IF
C                                             Have completed computation
C                                             for new iterate B.
 20   CALL SDNTP (H, 0, N, NQ, T, B, YH,  Y)
      FB = F(N, B, Y, IROOT)
      IF (FB .EQ. 0.E0) RETURN
      KOUNT = KOUNT + 1
C
C             Decide whether next step is interpolation or extrapolation
C
      IF (SIGN(1.0E0, FB) .EQ. SIGN(1.0E0, FC)) THEN
        C = A
        FC = FA
      END IF
      GO TO 10
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
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
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
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
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
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
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
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
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
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
C     Latest revision ---  19 MAR 1980
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
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
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
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
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
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
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
      subroutine caxpy ( n, ca, cx, incx, cy, incy )

c*********************************************************************72
c
cc CAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, integer N, the number of elements in CX and CY.
c
c    Input, complex CA, the multiplier of CX.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), the second vector.
c    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      complex ca
      complex cx(*)
      complex cy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( abs ( real ( ca ) ) + abs ( aimag ( ca ) ) .eq. 0.0 ) then
        return
      end if
c
c  Code for both increments equal to 1.
c
      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        do i = 1, n
          cy(i) = cy(i) + ca * cx(i)
        end do
c
c  Code for unequal increments or equal increments
c  not equal to 1.
c
      else

        if ( incx .lt. 0 ) then
          ix = ( - n + 1 ) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = ( - n + 1 ) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          cy(iy) = cy(iy) + ca * cx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      return
      end
      function cdotc ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CDOTC forms the conjugated  dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, complex CDOTC, the conjugated dot product of
c    the corresponding entries of CX and CY.
c
      implicit none

      complex cdotc
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      cdotc = ctemp
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
      end do

      cdotc = ctemp
      return
      end
      subroutine cgbfa(abd,lda,n,ml,mu,ipvt,info)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),info
      complex abd(lda,1)
c
c     cgbfa factors a complex band matrix by elimination.
c
c     cgbfa is usually called by cgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     complex(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that cgbsl will divide by zero if
c                     called.  use  rcond  in cgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,max0,min0,real
c
c     internal variables
c
      complex t
      integer i,icamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = (0.0e0,0.0e0)
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = (0.0e0,0.0e0)
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = icamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(abd(l,k)) .eq. 0.0e0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -(1.0e0,0.0e0)/abd(m,k)
            call cscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call caxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (cabs1(abd(m,n)) .eq. 0.0e0) info = n
      return
      end
      subroutine cgbsl(abd,lda,n,ml,mu,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),job
      complex abd(lda,1),b(1)
c
c     cgbsl solves the complex band system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cgbco or cgbfa.
c
c     on entry
c
c        abd     complex(lda, n)
c                the output from cgbco or cgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from cgbco or cgbfa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b , where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if cgbco has set rcond .gt. 0.0
c        or cgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call cgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg,min0
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call caxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call caxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = cdotc(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/conjg(abd(m,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + cdotc(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine cgefa(a,lda,n,ipvt,info)

c*********************************************************************72

      integer lda,n,ipvt(1),info
      complex a(lda,1)
c
c     cgefa factors a complex matrix by gaussian elimination.
c
c     cgefa is usually called by cgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cgeco) = (1 + 9/n)*(time for cgefa) .
c
c     on entry
c
c        a       complex(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that cgesl or cgedi will divide by zero
c                     if called.  use  rcond  in cgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,real
c
c     internal variables
c
      complex t
      integer icamax,j,k,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = icamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0e0,0.0e0)/a(k,k)
            call cscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call caxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0e0) info = n
      return
      end
      subroutine cgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      complex a(lda,1),b(1)
c
c     cgesl solves the complex system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by cgeco or cgefa.
c
c     on entry
c
c        a       complex(lda, n)
c                the output from cgeco or cgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from cgeco or cgefa.
c
c        b       complex(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if cgeco has set rcond .gt. 0.0
c        or cgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call cgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call cgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cdotc
c     fortran conjg
c
c     internal variables
c
      complex cdotc,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call caxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call caxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = cdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/conjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + cdotc(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine cscal ( n, ca, cx, incx )

c*********************************************************************72
c
cc CSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex CA, the multiplier.
c
c    Input/output, complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
      implicit none

      complex ca,cx(*)
      integer i,incx,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        cx(i) = ca*cx(i)
      end do
      return
c
c  code for increment equal to 1
c
   20 do i = 1,n
        cx(i) = ca*cx(i)
      end do

      return
      end
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    Input/output unit numbers.
c
c      I1MACH(1) = the standard input unit.
c      I1MACH(2) = the standard output unit.
c      I1MACH(3) = the standard punch unit.
c      I1MACH(4) = the standard error message unit.
c
c    Words.
c
c      I1MACH(5) = the number of bits per integer storage unit.
c      I1MACH(6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S digit base A form:
c
c      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A**S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit
c    base B form:
c
c      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
c
c    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
c
c      I1MACH(10) = B, the base.
c
c    Single precision
c
c      I1MACH(11) = T, the number of base B digits.
c      I1MACH(12) = EMIN, the smallest exponent E.
c      I1MACH(13) = EMAX, the largest exponent E.
c
c    Double precision
c
c      I1MACH(14) = T, the number of base B digits.
c      I1MACH(15) = EMIN, the smallest exponent E.
c      I1MACH(16) = EMAX, the largest exponent E.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 16.
c
c    Output, integer I1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      integer i1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      else if ( i == 1 ) then
        i1mach = 5
      else if ( i == 2 ) then
        i1mach = 6
      else if ( i == 3 ) then
        i1mach = 7
      else if ( i == 4 ) then
        i1mach = 6
      else if ( i == 5 ) then
        i1mach = 32
      else if ( i == 6 ) then
        i1mach = 4
      else if ( i == 7 ) then
        i1mach = 2
      else if ( i == 8 ) then
        i1mach = 31
      else if ( i == 9 ) then
        i1mach = 2147483647
      else if ( i == 10 ) then
        i1mach = 2
      else if ( i == 11 ) then
        i1mach = 24
      else if ( i == 12 ) then
        i1mach = -125
      else if ( i == 13 ) then
        i1mach = 128
      else if ( i == 14 ) then
        i1mach = 53
      else if ( i == 15 ) then
        i1mach = -1021
      else if ( i == 16 ) then
        i1mach = 1024
      else if ( 16 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      end if

      return
      end
      function icamax ( n, cx, incx )

c*********************************************************************72
c
cc ICAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, integer ICAMAX, the index of the element of maximum
c    absolute value.
c
      implicit none

      complex cx(*)
      integer icamax
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))

      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
      end do

      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do i = 2,n
        if ( smax .lt. cabs1(cx(i)) ) then
          icamax = i
          smax = cabs1(cx(i))
        end if
      end do

      return
      end
      function isamax ( n, sx, incx )

c*********************************************************************72
c
cc ISAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    25 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, real X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer ISAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      integer i
      integer incx
      integer isamax
      integer ix
      integer n
      real sx(*)
      real smax

      isamax = 0

      if ( n .lt. 1 .or. incx .le. 0 ) then
        return
      end if

      isamax = 1
      if ( n .eq. 1 ) then
        return
      end if

      if ( incx .ne. 1 ) then

        ix = 1
        smax = abs ( sx(1) )
        ix = ix + incx
        do i = 2, n
          if ( smax .lt. abs ( sx(ix) ) ) then
            isamax = i
            smax = abs ( sx(ix) )
          end if
          ix = ix + incx
        end do

      else

        smax = abs ( sx(1) )
        do i = 2, n
          if ( smax .lt. abs ( sx(i) ) ) then
            isamax = i
            smax = abs ( sx(i) )
          end if
        end do

      end if

      return
      end
      function r1mach ( i )

c*********************************************************************72
c
cc R1MACH returns single precision real machine constants.
c
c  Discussion:
c
c    Assume that single precision real numbers are stored with a mantissa
c    of T digits in base B, with an exponent whose value must lie
c    between EMIN and EMAX.  Then for values of I between 1 and 5,
c    R1MACH will return the following values:
c
c      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
c      R1MACH(3) = B**(-T), the smallest relative spacing.
c      R1MACH(4) = B**(1-T), the largest relative spacing.
c      R1MACH(5) = log10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 5.
c
c    Output, real R1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      real r1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      else if ( i == 1 ) then
        r1mach = 1.1754944E-38
      else if ( i == 2 ) then
        r1mach = 3.4028235E+38
      else if ( i == 3 ) then
        r1mach = 5.9604645E-08
      else if ( i == 4 ) then
        r1mach = 1.1920929E-07
      else if ( i == 5 ) then
        r1mach = 0.3010300E+00
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      end if

      return
      end
      subroutine saxpy ( n, sa, sx, incx, sy, incy )

c*********************************************************************72
c
cc SAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loop for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input, real X(*), the vector to be scaled and added to Y.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), the vector to which a multiple of X is to
c    be added.
c
c    Input, integer INCY, the increment between successive entries of Y.
c
      implicit none

      real sx(*),sy(*),sa
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sy(i) = sy(i) + sa*sx(i)
      end do
      if( n .lt. 4 ) return
   40 continue

      do i = m+1, n, 4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
      end do

      return
      end
      function scnrm2 ( n, x, incx )

c*********************************************************************72
c
cc SCNRM2 returns the euclidean norm of a complex vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
c            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Sven Hammarling
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SCNRM2, the norm of the vector.
c
      implicit none

      integer incx
      integer ix
      integer n
      real one
      parameter ( one = 1.0e+0 )
      real scnrm2
      complex x( * )
      real zero
      parameter ( zero = 0.0e+0 )


      real                  norm, scale, ssq, temp
      intrinsic             abs, aimag, real, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else
         scale = zero
         ssq   = one
c        The following loop is equivalent to this call to the LAPACK
c        auxiliary routine:
c        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( real( x( ix ) ).ne.zero )then
               temp = abs( real( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
            if( aimag( x( ix ) ).ne.zero )then
               temp = abs( aimag( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      scnrm2 = norm

      return
      end
      function sdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), one of the vectors to be multiplied.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input, real Y(*), one of the vectors to be multiplied.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Output, real SDOT, the dot product of X and Y.
c
      implicit none

      real sdot
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,n

      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      sdot = stemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + sx(i)*sy(i)
      end do

      if( n .lt. 5 ) go to 60
   40 continue

      do i = m+1, n, 5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
      end do

   60 sdot = stemp

      return
      end
      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),info
      real abd(lda,1)
c
c     sgbfa factors a real band matrix by elimination.
c
c     sgbfa is usually called by sgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     real(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgbsl will divide by zero if
c                     called.  use  rcond  in sgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c     fortran max0,min0
c
c     internal variables
c
      real t
      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0e0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0e0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = isamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0e0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0e0/abd(m,k)
            call sscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call saxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0e0) info = n
      return
      end
      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ml,mu,ipvt(1),job
      real abd(lda,1),b(1)
c
c     sgbsl solves the real band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgbco or sgbfa.
c
c     on entry
c
c        abd     real(lda, n)
c                the output from sgbco or sgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from sgbco or sgbfa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgbco has set rcond .gt. 0.0
c        or sgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c     fortran min0
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call saxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call saxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = sdot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine sgefa(a,lda,n,ipvt,info)

c*********************************************************************72

      integer lda,n,ipvt(1),info
      real a(lda,1)
c
c     sgefa factors a real matrix by gaussian elimination.
c
c     sgefa is usually called by sgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c
c     on entry
c
c        a       real(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that sgesl or sgedi will divide by zero
c                     if called.  use  rcond  in sgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sscal,isamax
c
c     internal variables
c
      real t
      integer isamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = isamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end
      subroutine sgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      real a(lda,1),b(1)
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call sgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas saxpy,sdot
c
c     internal variables
c
      real sdot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call saxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = sdot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      function snrm2 ( n, x, incx )

c*********************************************************************72
c
cc SNRM2 returns the euclidean norm of a real vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Sven Hammarling
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, real X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SNRM2, the Euclidean norm of X.
c
      implicit none

      integer                           incx, n
      real snrm2
      real                              x( * )

      real                  one         , zero
      parameter           ( one = 1.0e+0, zero = 0.0e+0 )

      integer               ix
      real                  absxi, norm, scale, ssq

      intrinsic             abs, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
c        The following loop is equivalent to this call to the LAPACK
c        auxiliary routine:
c        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
        end do
        norm  = scale * sqrt( ssq )
      end if

      snrm2 = norm
      return
      end
      subroutine sscal ( n, sa, sx, incx )

c*********************************************************************72
c
cc SSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input/output, real X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      real sa,sx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        sx(i) = sa*sx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sx(i) = sa*sx(i)
      end do
      if( n .lt. 5 ) return
   40 continue

      do i = m+1, n, 5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
      end do

      return
      end
