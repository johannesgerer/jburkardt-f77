C
C                          CS2TST
C                         11/20/98
C
C   This program computes interpolation errors using the
C scattered data package CSHEP2D for each of ten test
C functions and a 33 by 33 uniform grid of interpolation
C points in the unit square.
C
C   This program uses Subroutines TESTDT and TSTFN1 from
C ACM Algorithm SURVEY to generate a node set and and the
C test function values.
C
      INTEGER   NMAX, NRMAX, NI
      PARAMETER (NMAX=100, NRMAX=10, NI=33)
C
C Array storage:
C
      DOUBLE PRECISION X(NMAX), Y(NMAX), W(NMAX), RW(NMAX),
     .                 A(9,NMAX), P(NI), FT(NI,NI)
      INTEGER          LCELL(NRMAX,NRMAX), LNEXT(NMAX)
C
      DOUBLE PRECISION DEL, DUM, DX, DY, ERMAX, ERMEAN, PW,
     .                 RMAX, SSA, SSE, SSM, SUM, XMIN, YMIN
      DOUBLE PRECISION CS2VAL
      INTEGER          I, IER, J, K, KF, KFF, KFL, KS,
     .                 N, NC, NFUN, NP, NR, NSET, NW, NWMAX
C
C Data:
C
C NSET = Number of node sets.
C NFUN = Number of test functions.
C
      DATA  NSET/5/, NFUN/10/
C
C Input format:
C
  100 FORMAT (I2)
C
C Get a user-specified node set number KS.
C
    1 WRITE (*,110) NSET
  110 FORMAT (///13X,'CS2TST:  CSHEP2D Test Program'//
     .        5X,'Specify a node set number in the range 1',
     .           ' to ',I2,':'/)
      READ (*,100,END=999,ERR=1) KS
      IF (KS .LT. 1  .OR.  KS .GT. NSET) GO TO 1
C
C Copy N and the nodal coordinates for node set KS.
C
      CALL TESTDT (KS, N,X,Y)
      IF (N .LT. 10  .OR.  N .GT. NMAX) GO TO 20
C
C Allow the user to specify a range of function numbers.
C
      WRITE (*,120) NFUN
  120 FORMAT (//5X,'Specify the first test function ',
     .             '(1 to ',I2,'):'/)
      READ (*,100,ERR=1) KFF
      IF (KFF .LT. 1  .OR.  KFF .GT. NFUN) GO TO 1
      write ( *, * ) '  KFF = ', kff

      WRITE (*,130) KFF, NFUN
  130 FORMAT (//5X,'Specify the last test function (',
     .             I2,' to ',I2,'):'/)
      READ (*,100,ERR=1) KFL
      IF (KFL .LT. KFF  .OR.  KFL .GT. NFUN) GO TO 1
      write ( *, * ) '  KFL = ', kfl

      NFUN = KFL-KFF+1
C
C Input NC, NW, and NR from the console.
C
      NWMAX = MIN(40,N-1)
    2 WRITE (*,140) N
  140 FORMAT (//5X,'N =',I4//5X,
     .        'Specify the number of nodes NC for the ',
     .        'least squares fits.'/5X,'NC = 17 is ',
     .        'recommended.  NC GE 9.'/)
      READ (*,100,ERR=2) NC
      IF (NC .LT. 9  .OR.  NC .GT. NWMAX) GO TO 2
      write ( *, * ) '  NC = ', nc
C
    3 WRITE (*,150)
  150 FORMAT (///5X,'Specify the number of nodes NW for ',
     .        'the weights.  NW = 30 is'/5X,'recommended. ',
     .        ' 1  LE  NW  LE  MIN(40,N-1).')
      READ (*,100,ERR=2) NW
      IF (1 .GT. NW  .OR.  NW .GT. NWMAX) GO TO 3
      write ( *, * ) '  NW = ', nw
C
    4 WRITE (*,160) NRMAX
  160 FORMAT (///5X,'Specify the number of rows and column',
     .        's NR in the uniform grid'/5X,'of cells used',
     .        ' to locate nearest neighbors.  NR = Sqrt(N/',
     .        '3) is'/5X,'recommended.  1 LE NR LE ',I2)
      READ (*,100,ERR=3) NR
      IF (NR .LT. 1  .OR.  NR .GT. NRMAX) GO TO 4
      write ( *, * ) '  NR = ', nr
C
C Set up uniform grid points.
C
      DEL = 1./DBLE(NI-1)
      DO 5 I = 1,NI
        P(I) = DBLE(I-1)*DEL
    5   CONTINUE
C
C Initialize the average SSE/SSM value to zero.
C
      SSA = 0.
C
C Print a heading and loop on test functions.
C
      WRITE (*,200) KS, N, NI, NC, NW, NR
      DO 11 KF = KFF,KFL
C
C   Compute true function values at the nodes.
C
        DO 6 K = 1,N
          CALL TSTFN1 (KF,X(K),Y(K),0, W(K),DUM,DUM)
    6     CONTINUE
C
C   Compute true function values FT on the uniform grid, and
C     accumulate the sum of values SUM and sum of squared
C     values SSM.
C
        SUM = 0.
        SSM = 0.
        DO 8 I = 1,NI
          DO 7 J = 1,NI
            CALL TSTFN1 (KF,P(I),P(J),0, FT(I,J),DUM,DUM)
            SUM = SUM + FT(I,J)
            SSM = SSM + FT(I,J)**2
    7       CONTINUE
    8     CONTINUE
C
C   Compute the sum of squared deviations from the mean SSM.
C
        SSM = SSM - SUM*SUM/DBLE(NI*NI)
C
C   Compute parameters A and RW defining the interpolant.
C
        CALL CSHEP2 (N,X,Y,W,NC,NW,NR, LCELL,LNEXT,XMIN,
     .               YMIN,DX,DY,RMAX,RW,A,IER)
        IF (IER .NE. 0) GO TO 21
C
C   Compute interpolation errors.
C
        ERMEAN = 0.
        ERMAX = 0.
        SSE = 0.
        DO 10 I = 1,NI
          DO 9 J = 1,NI
            PW = CS2VAL (P(I),P(J),N,X,Y,W,NR,LCELL,LNEXT,
     .                   XMIN,YMIN,DX,DY,RMAX,RW,A) -
     .           FT(I,J)
            ERMEAN = ERMEAN + ABS(PW)
            ERMAX = MAX(ERMAX,ABS(PW))
            SSE = SSE + PW*PW
    9       CONTINUE
   10     CONTINUE
        NP = NI*NI
        ERMEAN = ERMEAN/DBLE(NP)
        SSE = SSE/SSM
        SSA = SSA + SSE
        WRITE (*,210) KF, ERMAX, ERMEAN, SSE
   11   CONTINUE
C
C Print the average SSE/SSM value (averaged over the test
C   functions).
C
      WRITE (*,220) SSA/DBLE(NFUN)
      go to 1
C
C N is outside its valid range.
C
   20 WRITE (*,500) N, NMAX
      STOP
C
C Error in CSHEP2.
C
   21 IF (IER .EQ. 2) WRITE (*,510)
      IF (IER .EQ. 3) WRITE (*,520)
      STOP

999   continue
      write ( *, * ) ' '
      write ( *, '(a)' ) 'TOMS792_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      stop
C
C Print formats:
C
  200 FORMAT (30X,'CS2TST Output:'//
     .        1X,16X,'CSHEP2D Interpolation Errors for ',
     .        'N nodes and an'/
     .        1X,6X,'NI by NI Uniform Grid of Interpolation',
     .        'ion Points in the Unit Square'//1X,
     .        6X,'Node set ',I2,4X,'N =',I4,4X,'NI = ',I2,
     .        4X,'NC = ',I2,4X,'NW = ',I2,4X,'NR = ',I2///
     .        1X,16X,'Function',4X,'Max Error',4X,
     .        'Mean Error',4X,'SSE/SSM'/)
  210 FORMAT (1X,19X,I2,9X,F7.4,6X,F8.5,2X,F9.6)
  220 FORMAT (//1X,11X,'Average SSE/SSM over the test ',
     .        'functions = ',F9.6)
C
C Error message formats:
C
  500 FORMAT (///1X,10X,'*** Error in data -- N = ',I4,
     .        ', Maximum value =',I4,' ***')
  510 FORMAT (///1X,14X,'*** Error in CSHEP2 -- duplicate ',
     .        'nodes encountered ***')
  520 FORMAT (///1X,14X,'*** Error in CSHEP2 -- all nodes ',
     .        'are collinear ***')
      END
