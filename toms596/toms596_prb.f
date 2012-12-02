      program main

C***********************************************************************
C
C  TEST MAIN PROGRAM FOR CONTINUATION CODE
C  AIRCRAFT STABILITY PROBLEM
C
C  RUN WITH X(6) = VAL1 =  -.05, -.008, 0.0, 0.1
C
C  SEEKING LIMIT POINTS IN X(7).
C
C  X(1) = ROLL
C  X(2) = PITCH
C  X(3) = YAW
C  X(4) = ANGLE OF ATTACK
C  X(5) = SIDESLIP
C  X(6) = ELEVATOR
C  X(7) = AILERON
C  X(8) = RUDDER
C
      EXTERNAL FXAIR, FPAIR, FSOLVE
      REAL RWORK(104)
      INTEGER IPIVOT(8)
      COMMON /COEF/ BARRAY(5,8)
      COMMON /CONTRL/ IFIX1, IFIX2, VAL1, VAL2
      COMMON /COUNT1/ ICRSL, ITNSL, NSTCR, NCNCR, NTRCR, NLMCR, NLMRT
      COMMON /COUNT2/ IFEVAL, IPEVAL, ISOLVE, NREDXF, NRDSUM, NCORXR,
     * NCRSUM
      COMMON /DETEXP/ IEXTXL, IEXTXC, IEXTXR, IEXNXR
      COMMON /DETMAN/ DETTXL, DETTXC, DETTXR, DETNXR
      COMMON /OUTPUT/ IWRITE                                            MAN  310
      COMMON /POINT/ CURVCF, CORDXF, ALPHLC, HSECLC, FNRM               MAN  320
      IWRITE = 1                                                        MAN  330
      LUNIT = 6                                                         MAN  340
      MAXSTP = 25                                                       MAN  350
C                                                                       MAN  360
C  IFIX1 AND IFIX2 RECORD THE INDICES OF X WHICH ARE TO BE HELD FIXED.  MAN  370
C  VAL1 AND VAL2 ARE THE VALUES AT WHICH THESE QUANTITIES ARE HELD.     MAN  380
C                                                                       MAN  390
      IFIX1 = 6                                                         MAN  400
      IFIX2 = 8                                                         MAN  410
      VAL1 = -.05                                                       MAN  420
      VAL2 = 0.0                                                        MAN  430
C                                                                       MAN  440
C  SET PITCON ARGUMENTS                                                 MAN  450
C                                                                       MAN  460
      NVAR = 8                                                          MAN  470
      LIM = 7                                                           MAN  480
      IT = 0                                                            MAN  490
      XIT = 0.0                                                         MAN  500
      KSTEP2 = -3                                                       MAN  510
      KSTEP1 = -2                                                       MAN  520
      KSTEP = -1                                                        MAN  530
      IPC = 7                                                           MAN  540
      IPCFIX = 0                                                        MAN  550
      DIRIPC = -1.0                                                     MAN  560
      HTANCF = .3                                                       MAN  570
      IRET = 0                                                          MAN  580
      MODCON = 0                                                        MAN  590
      HMAX = 100.0                                                      MAN  600
      HMIN = .001                                                       MAN  610
      HFACT = 3.0                                                       MAN  620
      ABSERR = .00001                                                   MAN  630
      RELERR = .00001                                                   MAN  640
      RWORK(1) = 0.0                                                    MAN  650
      RWORK(2) = 0.0                                                    MAN  660
      RWORK(3) = 0.0                                                    MAN  670
      RWORK(4) = 0.0                                                    MAN  680
      RWORK(5) = 0.0                                                    MAN  690
      RWORK(6) = 0.0                                                    MAN  700
      RWORK(7) = 0.0
      RWORK(8) = 0.0
      RWORK(IFIX1) = VAL1
      RWORK(IFIX2) = VAL2
      ISIZE = 104
      NROW = NVAR
      NCOL = NVAR
C
C  SET LOCATION POINTERS
C
      JXR = 0
      IXR = 1
      JXC = JXR + NVAR
      IXC = IXR + NVAR
      JXF = JXC + NVAR
      IXF = IXC + NVAR
      JTL = JXF + NVAR
      ITL = IXF + NVAR
      JTC = JTL + NVAR
      ITC = ITL + NVAR
      JFP = JTC + NVAR
      IFP = ITC + NVAR
C
C  SET B MATRIX
C
      CALL STORE(BARRAY)
      WRITE (LUNIT,99987) RWORK(6), RWORK(7), RWORK(8)
      WRITE (LUNIT,99988) ABSERR, RELERR, HTANCF, HFACT, HMIN, HMAX,
     * MODCON
      IF (LIM.NE.0) WRITE (LUNIT,99984) LIM
      IF (IPCFIX.NE.0) WRITE (LUNIT,99983) IPC                          MAN 1010
C                                                                       MAN 1020
C  BEGIN TRACING OUT THE CURVE                                          MAN 1030
C                                                                       MAN 1040
      DO 20 I=1,MAXSTP                                                  MAN 1050
        CALL PITCON(NVAR, LIM, IT, XIT, KSTEP, IPC, IPCFIX, DIRIPC,     MAN 1060
     *   HTANCF, IRET, MODCON, IPIVOT, HMAX, HMIN, HFACT, ABSERR,       MAN 1070
     *   RELERR, RWORK, ISIZE, NROW, NCOL, FXAIR, FPAIR, FSOLVE, LUNIT) MAN 1080
        IF (IRET.LT.(-1)) GO TO 30                                      MAN 1090
        IF (IRET.LT.0) GO TO 20                                         MAN 1100
        IF (IRET.EQ.1) GO TO 40                                         MAN 1110
        IF (IRET.EQ.2) GO TO 10                                         MAN 1120
C                                                                       MAN 1130
C  NORMAL CONTINUATION POINT RETURN                                     MAN 1140
C                                                                       MAN 1150
        KSTEP1 = KSTEP - 1                                              MAN 1160
        KSTEP2 = KSTEP1 - 1                                             MAN 1170
        IF (KSTEP.GT.1) WRITE (LUNIT,99996) KSTEP2, KSTEP1, HSECLC      MAN 1180
        IHI = JTC + NVAR                                                MAN 1190
        WRITE (LUNIT,99995) KSTEP1                                      MAN 1200
        WRITE (LUNIT,99982) (RWORK(J),J=ITC,JFP)                        MAN 1210
        WRITE (LUNIT,99994) KSTEP, HTANCF                               MAN 1220
        WRITE (LUNIT,99993) KSTEP                                       MAN 1230
        WRITE (LUNIT,99982) (RWORK(J),J=IXR,JXC)                        MAN 1240
        WRITE (LUNIT,99992) FNRM                                        MAN 1250
        GO TO 20                                                        MAN 1260
C                                                                       MAN 1270
C  LIMIT POINT FOUND, PRINT OUT AND CONTINUE                            MAN 1280
C                                                                       MAN 1290
   10   WRITE (LUNIT,99991)
        WRITE (LUNIT,99986)
        WRITE (LUNIT,99982) (RWORK(J),J=IXR,JXC)
        WRITE (LUNIT,99989)
        WRITE (LUNIT,99982) (RWORK(J),J=ITL,JTC)
        WRITE (LUNIT,99992) FNRM
   20 CONTINUE
C
C  LOOP COMPLETED WITHOUT REACHING TARGET POINT
C
      WRITE (LUNIT,99999)
      GO TO 50
C
C  ERROR RETURN FROM CONTINUATION CODE
C
   30 WRITE (LUNIT,99998) IRET
      GO TO 50
C
C  TARGET POINT REACHED
C
   40 WRITE (LUNIT,99997)
      WRITE (LUNIT,99986)
      WRITE (LUNIT,99982) (RWORK(J),J=IXR,JXC)
      WRITE (LUNIT,99992) FNRM
   50 WRITE (LUNIT,99990) NROW, NCOL, KSTEP, NCRSUM, NRDSUM, IFEVAL,
     * IPEVAL, ISOLVE
      WRITE (LUNIT,99985) ICRSL, ITNSL, NSTCR, NCNCR, NTRCR, NLMCR,
     * NLMRT
      STOP
99999 FORMAT (47H0PROGRAM TOOK FULL NUMBER OF CONTINUATION STEPS)
99998 FORMAT (24H0ERROR RETURN WITH IRET=, I4)
99997 FORMAT (19H0TARGET FOUND, STOP)
99996 FORMAT (17H0SECANT STEP FROM, I3, 4H TO , I3, 3H = , G14.7)
99995 FORMAT (18H TANGENT AT POINT , I3, 5H WAS )
99994 FORMAT (16H TO REACH POINT , I3, 22H TANGENT STEPSIZE WAS , G14.7)
99993 FORMAT (7H POINT , I3, 4H IS )
99992 FORMAT (15H FUNCTION NORM=, G14.7)
99991 FORMAT (18H0LIMIT POINT FOUND)
99990 FORMAT (17H0MATRIX SIZE     , I4, 1HX, I4/17H STEPS TAKEN     ,
     * I4/17H NEWTON STEPS    , I4/17H STEP REDUCTIONS , I4/9H FUNCTION,
     * 8H CALLS  , I4/17H PRIME CALLS     , I4/17H SOLVE CALLS     , I4)
99989 FORMAT (9H TANGENT )
99988 FORMAT (43H0UNIVERSITY OF PITTSBURGH CONTINUATION CODE/8H ABSERR=,
     * G14.7, 8H RELERR=, G14.7/8H H=     , G14.7, 8H HFACT= ,
     * G14.7/8H HMIN=  , G14.7, 8H HMAX=  , G14.7/18H NEWTON METHOD OPT,
     * 4HION=, I2)
99987 FORMAT (27H1AIRCRAFT STABILITY PROBLEM/10H ELEVATOR=,
     * G14.7/10H AILERON =, G14.7/10H RUDDER  =, G14.7)
99986 FORMAT (9H POINT   )
99985 FORMAT (8H0CALLS -/40H SOLVE CALLED FOR CORRECTION            ,
     * I3/40H SOLVE CALLED FOR TANGENT               , I3/10H0CORRECTOR,
     * 30H CALLED FOR STARTING POINT    , I3/23H CORRECTOR CALLED FOR C,
     * 17HONTINUATION POINT, I3/36H CORRECTOR CALLED FOR TARGET POINT
     * 4H    , I3/40H CORRECTOR CALLED FOR LIMIT POINT       ,
     * I3/40H0ROOTFINDER CALLED FOR LIMIT POINT      , I3)
99984 FORMAT (33H0LIMIT POINTS SOUGHT FOR VARIABLE, I3)
99983 FORMAT (36H CONTINUATION VARIABLE FORCED TO BE , I3)
99982 FORMAT (1X, 5G14.7)
      END
      SUBROUTINE STORE(BARRAY)

C***********************************************************************
C
      REAL BARRAY(5,8)
      BARRAY(1,1) = -3.933
      BARRAY(2,1) = 0.0
      BARRAY(3,1) = 0.002
      BARRAY(4,1) = 0.0
      BARRAY(5,1) = 0.0
      BARRAY(1,2) = 0.107
      BARRAY(2,2) = -0.987
      BARRAY(3,2) = 0.0
      BARRAY(4,2) = 1.0
      BARRAY(5,2) = 0.0
      BARRAY(1,3) = 0.126
      BARRAY(2,3) = 0.0
      BARRAY(3,3) = -0.235
      BARRAY(4,3) = 0.0
      BARRAY(5,3) = -1.0
      BARRAY(1,4) = 0.0
      BARRAY(2,4) = -22.95
      BARRAY(3,4) = 0.0
      BARRAY(4,4) = -1.0
      BARRAY(5,4) = 0.0
      BARRAY(1,5) = -9.99
      BARRAY(2,5) = 0.0
      BARRAY(3,5) = 5.67
      BARRAY(4,5) = 0.0
      BARRAY(5,5) = -0.196
      BARRAY(1,6) = 0.0
      BARRAY(2,6) = -28.37
      BARRAY(3,6) = 0.0
      BARRAY(4,6) = -0.168
      BARRAY(5,6) = 0.0
      BARRAY(1,7) = -45.83
      BARRAY(2,7) = 0.0
      BARRAY(3,7) = -0.921
      BARRAY(4,7) = 0.0
      BARRAY(5,7) = -0.0071
      BARRAY(1,8) = -7.64
      BARRAY(2,8) = 0.0
      BARRAY(3,8) = -6.51
      BARRAY(4,8) = 0.0
      BARRAY(5,8) = 0.0

      RETURN
      END
      SUBROUTINE FXAIR(NVAR, X, FX)
C
C***********************************************************************
C
C  THE FUNCTION IS OF THE FOLLOWING FORM
C
C  FOR INDICES I=1 THROUGH 5,
C
C  F(I)=SUM ON J AND K OF  B(I,J)*X(J)+PHI(I,J,K)*X(J)*X(K)
C
C  F(6)=X(IFIX1)-VAL1
C  F(7)=X(IFIX2)-VAL2
C
      REAL X(NVAR), FX(NVAR)
      COMMON /COEF/ BARRAY(5,8)
      COMMON /CONTRL/ IFIX1, IFIX2, VAL1, VAL2
C
C  COMPUTE LINEAR TERMS
C
      DO 20 I=1,5
        FX(I) = 0.0
        DO 10 J=1,8
          FX(I) = FX(I) + BARRAY(I,J)*X(J)
   10   CONTINUE
   20 CONTINUE
C
C  COMPUTE NONLINEAR TERMS
C
      PHI = -.727*X(2)*X(3) + 8.39*X(3)*X(4) - 684.4*X(4)*X(5) +
     * 63.5*X(4)*X(7)
      FX(1) = FX(1) + PHI
      PHI = .949*X(1)*X(3) + .173*X(1)*X(5)
      FX(2) = FX(2) + PHI
      PHI = -.716*X(1)*X(2) - 1.578*X(1)*X(4) + 1.132*X(4)*X(7)
      FX(3) = FX(3) + PHI
      PHI = -X(1)*X(5)
      FX(4) = FX(4) + PHI
      PHI = X(1)*X(4)
      FX(5) = FX(5) + PHI
C
C  SET FUNCTION VALUES FOR TWO FIXED VARIABLES
C
      FX(6) = X(IFIX1) - VAL1
      FX(7) = X(IFIX2) - VAL2

      RETURN
      END
      SUBROUTINE FPAIR(NVAR, X, FPRYM, NROW, NCOL)

C***********************************************************************
C
      REAL X(NVAR), FPRYM(NROW,NCOL)
      COMMON /COEF/ BARRAY(5,8)
      COMMON /CONTRL/ IFIX1, IFIX2, VAL1, VAL2
C
C  COMPUTE TERMS FROM LINEAR PART OF FUNCTION
C
      DO 20 I=1,5
        DO 10 J=1,8
          FPRYM(I,J) = BARRAY(I,J)
   10   CONTINUE
   20 CONTINUE
      DO 40 I=6,7
        DO 30 J=1,8
          FPRYM(I,J) = 0.0
   30   CONTINUE
   40 CONTINUE
      FPRYM(6,IFIX1) = 1.0
      FPRYM(7,IFIX2) = 1.0
C
C  COMPUTE TERMS FROM NONLINEAR PART OF FUNCTION
C
      FPRYM(1,2) = FPRYM(1,2) - .727*X(3)
      FPRYM(1,3) = FPRYM(1,3) - .727*X(2) + 8.39*X(4)
      FPRYM(1,4) = FPRYM(1,4) + 8.39*X(3) - 684.4*X(5) + 63.5*X(7)
      FPRYM(1,5) = FPRYM(1,5) - 684.4*X(4)
      FPRYM(1,7) = FPRYM(1,7) + 63.5*X(4)
      FPRYM(2,1) = FPRYM(2,1) + .949*X(3) + .173*X(5)
      FPRYM(2,3) = FPRYM(2,3) + .949*X(1)
      FPRYM(2,5) = FPRYM(2,5) + .173*X(1)
      FPRYM(3,1) = FPRYM(3,1) - .716*X(2) - 1.578*X(4)
      FPRYM(3,2) = FPRYM(3,2) - .716*X(1)
      FPRYM(3,4) = FPRYM(3,4) - 1.578*X(1) + 1.132*X(7)
      FPRYM(3,7) = FPRYM(3,7) + 1.132*X(4)
      FPRYM(4,1) = FPRYM(4,1) - X(5)
      FPRYM(4,5) = FPRYM(4,5) - X(1)
      FPRYM(5,1) = FPRYM(5,1) + X(4)
      FPRYM(5,4) = FPRYM(5,4) + X(1)

      RETURN
      END
