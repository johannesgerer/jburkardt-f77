*
      SUBROUTINE BICG( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                 MATVECTRANS, PSOLVE, PSOLVETRANS, INFO )
*
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
*
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS
*     ..
*
*  Purpose
*  =======
*
*  BiCG solves the linear system Ax = b using the
*  BiConjugate Gradient iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X      (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,6).
*          Workspace for residual, direction vector, etc.
*          Note that Z and Q, and ZTLD and QTLD share workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  MATVECTRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A'*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A' is the tranpose of a matrix A. Vector x must remain 
*          unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVECTRANS( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  PSOLVETRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a 
*          matrix M. Vector b must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANS( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  ==============================================================
*
*
*This variable used to communicate requests between BiCG() and 
* BiCG_CGREVCOM()
*BiCG -> BiCGREVCOM: 1 = init, 
*                    2 = use saved state to resume flow.
*BiCGREVCOM -> BiCG: -1 = done, return to main, 
*                     1 = matvec using SCLR1/2, NDX1/2 
*                     2 = transpose-matvec using SCLR1/2, NDX1/2 
*                     3 = solve using NDX1/2
*                     4 = transpose-solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from BiCGREVCOM to BiCG.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         BICGREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == RTLD; 3 == Z; 4 == ZTLD; 5 == P; 6 == PTLD; 
*     7 == Q; 8 == QTLD; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL BICGREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from REVCOM we use the table (BiCGREVCOM -> BiCG)
*     to decode IJOB.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .eq. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 2) THEN
*        call transpose-matvec.
         CALL MATVECTRANS(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 3) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 4) THEN
*        call transpose-solve.
         CALL PSOLVETRANS(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 5) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 6) THEN
*        do stopping test 2
*        if FirstTime, set info  so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
         FTFLG = .FALSE.
      ENDIF
*
*     Done what revcom asked, set IJOB & go back to it.
*
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of BICG
*
      END
*
      SUBROUTINE BICGREVCOM( N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                       NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( LDW,* )
*
*     ..
*
*  Purpose
*  =======
*
*  BiCG solves the linear system Ax = b using the
*  BiConjugate Gradient iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X      (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,6).
*          Workspace for residual, direction vector, etc.
*          Note that Z and Q, and ZTLD and QTLD share workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become
*                                       orthogonal.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SDOT, SNRM2,
*  ==============================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, Z, ZTLD, P, PTLD, Q, QTLD, MAXIT,
     $                   NEED1, NEED2
      REAL               TOL, ALPHA, BETA, RHO, RHO1, BNRM2, RHOTOL, 
     $                   GETBREAK, SDOT, SNRM2
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     ..
*     .. External Routines ..
      EXTERNAL           SAXPY, SCOPY, SDOT, SNRM2
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
* init.
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      Z    = 3
      ZTLD = 4
      P    = 5
      PTLD = 6
      Q    = 3
      QTLD = 4
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((ZTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((QTLD - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((ZTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((QTLD - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set breakdown parameters.
*
      RHOTOL = GETBREAK()
* 
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ZERO, WORK(1,R) )
*        using WORK[RTLD] as temp
*********CALL SCOPY( N, X, 1, WORK(1,RTLD), 1 )
         SCLR1 = -ONE
         SCLR2 = ZERO
         NDX1 = ((RTLD - 1) * LDW) + 1
         NDX2 = ((R    - 1) * LDW) + 1
         RLBL = 2
         IJOB = 5
         RETURN
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
*
      ENDIF
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient iteration.
*
         ITER = ITER + 1
*
*        Compute direction vectors PK and PTLD.
*
*********CALL PSOLVE( WORK(1,Z), WORK(1,R) )
*
         NDX1 = ((Z - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
         RLBL = 3
         IJOB = 3
         RETURN
*****************
 3       CONTINUE
*****************
*
*********CALL PSOLVETRANS( WORK(1,ZTLD), WORK(1,RTLD) )
*
         NDX1 = ((ZTLD - 1) * LDW) + 1
         NDX2 = ((RTLD - 1) * LDW) + 1
         RLBL = 4
         IJOB = 4
         RETURN
*****************
 4       CONTINUE
*****************
*
*
         RHO = SDOT( N, WORK(1,Z), 1, WORK(1,RTLD), 1 )
         IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*
         IF ( ITER.GT.1 ) THEN
            BETA = RHO / RHO1
            CALL SAXPY( N, BETA, WORK(1,P), 1, WORK(1,Z), 1 )
            CALL SAXPY( N, BETA, WORK(1,PTLD), 1, WORK(1,ZTLD), 1 )
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
            CALL SCOPY( N, WORK(1,ZTLD), 1, WORK(1,PTLD), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
            CALL SCOPY( N, WORK(1,ZTLD), 1, WORK(1,PTLD), 1 )
         ENDIF
*
*********CALL MATVEC( ONE, WORK(1,P), ZERO, WORK(1,Q) )
*
         SCLR1 = ONE
         SCLR2 = ZERO
         NDX1 = ((P - 1) * LDW) + 1
         NDX2 = ((Q - 1) * LDW) + 1
         RLBL = 5
         IJOB = 1
         RETURN
*****************
 5       CONTINUE
*****************
*
*********CALL MATVECTRANS( ONE, WORK(1,PTLD), ZERO, WORK(1,QTLD) )
*
         SCLR1 = ONE
         SCLR2 = ZERO
         NDX1 = ((PTLD - 1) * LDW) + 1
         NDX2 = ((QTLD - 1) * LDW) + 1
         RLBL = 6
         IJOB = 2
         RETURN
*****************
 6       CONTINUE
*****************
         ALPHA = RHO / SDOT( N, WORK(1,PTLD), 1, WORK(1,Q), 1 )
*
*        Compute current solution vector x.
*
         CALL SAXPY( N, ALPHA, WORK(1,P), 1, X, 1 )
*
*        Compute residual vector rk, find norm,
*        then check for tolerance.
*
         CALL SAXPY( N, -ALPHA, WORK(1,Q), 1, WORK(1,R), 1 )
*
*********RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID.LE.TOL  ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 7
         IJOB = 6
         RETURN
*
*****************
 7       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         CALL SAXPY( N, -ALPHA, WORK(1,QTLD), 1, WORK(1,RTLD), 1 )
         RHO1 = RHO
*
      GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      INFO = -10
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGREVCOM
*
      END





*
      SUBROUTINE BICGSTAB( N, B, X, WORK, LDW, ITER, RESID, MATVEC,
     $                     PSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector b.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  ==============================================================
*     .. Local Scalars ..
*
*This variable used to communicate requests between BiCGSTAB() 
*  and BiCGSTABREVCOM()
*BiCGSTAB -> BiCGSTABREVCOM: 1 = init, 
*                            2 = use saved state to resume flow.
*BiCGSTABREVCOM -> BiCGSTAB: -1 = done, return to main, 
*                             1 = matvec using SCLR1/2, NDX1/2 
*                             2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from BiCGSTABREVCOM to BiCGSTAB.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         BICSTABGREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == RTLD; 3 == P; 4 == V; 5 == T; 6 == PHAT; 7 == SHAT;
*     8 == S; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL BICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                    NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from REVCOM() we use the table (REVCOM -> BiCGSTAB)
*     to decode IJOB.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .eq. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 2) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 3) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 4) THEN
*        do stopping test 2
*        if first time, set INFO so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
         FTFLG = .FALSE.
      ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
          GOTO 1
*
*     come here to terminate
*
 2    CONTINUE
*
      RETURN
*
*     End of BICGSTAB
*
      END
*
      SUBROUTINE BICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                    NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*  ==============================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      REAL               TOL, ALPHA, BETA, RHO, RHO1, BNRM2, OMEGA, 
     $                   RHOTOL, OMEGATOL, GETBREAK, SDOT, SNRM2
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Functions ..
      EXTERNAL           GETBREAK, SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = GETBREAK()
      OMEGATOL = GETBREAK()
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL SCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
      ENDIF
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2 .EQ. ZERO ) BNRM2 = ONE
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
         ITER = ITER + 1
*
         RHO = SDOT( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
         IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*
*        Compute vector P.
*
         IF ( ITER.GT.1 ) THEN
            BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
            CALL SAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
            CALL SSCAL( N, BETA, WORK(1,P), 1 )
            CALL SAXPY( N, ONE, WORK(1,R), 1, WORK(1,P), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
         ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((P    - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
         ALPHA = RHO / SDOT( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
*
*        Early check for tolerance.
*
         CALL SAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
         CALL SCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
         IF ( SNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
            CALL SAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
            RESID = SNRM2( N, WORK(1,S), 1 ) / BNRM2
            GO TO 30
         ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
            NDX1 = ((SHAT - 1) * LDW) + 1
            NDX2 = ((S    - 1) * LDW) + 1
*           Prepare for return & return
            RLBL = 5
            IJOB = 2
            RETURN
*
*****************
 5          CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
            NDX1 = ((SHAT - 1) * LDW) + 1
            NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
            SCLR1 = ONE
            SCLR2 = ZERO
            RLBL = 6
            IJOB = 1
            RETURN
*
*****************
 6          CONTINUE
*****************
*
            OMEGA = SDOT( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $              SDOT( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
            CALL SAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
            CALL SAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*
*           Compute residual R, check for tolerance.
*
            CALL SAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
            NDX1 = NEED1
            NDX2 = NEED2
*           Prepare for resumption & return
            RLBL = 7
            IJOB = 4
            RETURN
*
*****************
 7          CONTINUE
*****************
            IF( INFO.EQ.1 ) GO TO 30
*
            IF ( ITER.EQ.MAXIT ) THEN
               INFO = 1
               GO TO 20
            ENDIF
*
         ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*
      SUBROUTINE CG( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $               PSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ),  WORK( * )
*     ..
*     .. Subroutine Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  CG solves the linear system Ax = b using the
*  Conjugate Gradient iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*  --Done in CGREVCOM.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*          Set by CGREVCOM.
*
*  WORK    (workspace) REAL             array, dimension ( * ).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*          Set by CGREVCOM.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*          Set by CGREVCOM.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*          Set by CGREVCOM()
*  ============================================================
*
*     ..
*     .. Local Scalars ..
*This variable used to communicate requests between CG() and CGREVCOM()
*CG -> CGREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*CGREVCOM -> CG: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from CGREVCOM to CG.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         CGREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == Z; 3 == P; 4 == Q; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL CGREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $              NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from CGREVCOM() we use the table (CGREVCOM -> CG)
*     to figure out what is reqd.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .EQ. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 2) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .EQ. 3) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 4) THEN
*        do stopping test 2
*        if first time, set INFO so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
         FTFLG = .FALSE.
      ENDIF
*
*     done what revcom asked, set IJOB & go back to it.
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE

      RETURN
*
*     End of CG
*
      END
*
      SUBROUTINE CGREVCOM( N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ),  WORK( LDW,* )
*
*     (output) for matvec and solve. These index into WORK[]
*     INTEGER NDX1, NDX2
*     ..
*
*  Purpose
*  =======
*
*  CG solves the linear system Ax = b using the
*  Conjugate Gradient iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,4).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter.
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SDOT, SNRM2
*  ============================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER          MAXIT, R, Z, P, Q, NEED1, NEED2
      REAL             ALPHA, BETA, RHO, RHO1, SDOT, SNRM2,
     $                 TOL
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Routines ..
      EXTERNAL         SAXPY, SCOPY, SDOT, SNRM2
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
* init.
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL = RESID
*
*     Alias workspace columns.
*
      R = 1
      Z = 2
      P = 3
      Q = 4
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN

*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*
*        Set args for revcom return
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 30
      ENDIF
*
      ITER = 0
*
   10 CONTINUE
*
*        Perform Preconditioned Conjugate Gradient iteration.
*
         ITER = ITER + 1
*
*        Preconditioner Solve.
*
*********CALL PSOLVE( WORK(1,Z), WORK(1,R) )
*
         NDX1 = ((Z - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         RHO = SDOT( N, WORK(1,R), 1, WORK(1,Z), 1 )
*
*        Compute direction vector P.
*
         IF ( ITER.GT.1 ) THEN
            BETA = RHO / RHO1
            CALL SAXPY( N, BETA, WORK(1,P), 1, WORK(1,Z), 1 )
*
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
         ENDIF
*
*        Compute scalar ALPHA (save A*P to Q).
*
*********CALL MATVEC( ONE, WORK(1,P), ZERO, WORK(1,Q) )
*
         NDX1 = ((P - 1) * LDW) + 1
         NDX2 = ((Q - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
         ALPHA =  RHO / SDOT( N, WORK(1,P), 1, WORK(1,Q), 1 )
*
*        Compute current solution vector X.
*
         CALL SAXPY( N, ALPHA, WORK(1,P), 1, X, 1 )
*
*        Compute residual vector R, find norm,
*        then check for tolerance.
*
         CALL SAXPY( N, -ALPHA,  WORK(1,Q), 1, WORK(1,R), 1 )
*
*********RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 5
         IJOB = 4
         RETURN
*
*****************
 5       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         RHO1 = RHO
*
         GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of CGREVCOM
*
      END
*
      SUBROUTINE CGS( N, B, X, WORK, LDW, ITER, RESID, MATVEC, PSOLVE,
     $                INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  CGS solves the linear system Ax = b using the
*  Conjugate Gradient Squared iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector. The user should be warned that for
*          this particular algorithm, an initial guess close to
*          the actual solution can result in divergence.
*          On exit, the iterated solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors PHAT and QHAT, and UHAT and VHAT share 
*          the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On ouput, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product 
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged. 
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit.
*          >  0: Convergence not achieved. This will be set
*                to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occured
*                during iteration.
*
*                Illegal parameter: 
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If RHO become smaller than some tolerance,
*                   the program will terminate. Here we check 
*                   against tolerance BREAKTOL.
*
*                   -10: RHO < BREAKTOL: RHO and RTLD have become
*                                        orthogonal.
*
*  =============================================================
*     .. Local Scalars ..
*
*This variable used to communicate requests between CGS() and CGSREVCOM()
*CGS -> CGSREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*CGSREVCOM -> CGS: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from CGSREVCOM to CGS.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         CGSREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == RTLD; 3 == P; 4 == PHAT; 5 == Q; 6 == QHAT;
*     7 == U; 8 == UHAT; 9 == VHAT; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL CGSREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to figure out what is reqd.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call solve.
            CALL PSOLVE(WORK(NDX1), WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 4) THEN
*           do stopping test 2
*           if FirstTime, set info  so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
            FTFLG = .FALSE.
         ENDIF
*
*        done what revcom asked, set IJOB & go back to it.
         IJOB = 2
         GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of CGS
*
      END






*
      SUBROUTINE CGSREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  CGS solves the linear system Ax = b using the
*  Conjugate Gradient Squared iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector. The user should be warned that for
*          this particular algorithm, an initial guess close to
*          the actual solution can result in divergence.
*          On exit, the iterated solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors PHAT and QHAT, and UHAT and VHAT share 
*          the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On ouput, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit.
*          >  0: Convergence not achieved. This will be set
*                to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occured
*                during iteration.
*
*                Illegal parameter: 
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If RHO become smaller than some tolerance,
*                   the program will terminate. Here we check 
*                   against tolerance BREAKTOL.
*
*                   -10: RHO < BREAKTOL: RHO and RTLD have become
*                                        orthogonal.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:    SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*  =============================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER        ( ONE = 1.0E+0 , ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, Q, QHAT, U, UHAT, VHAT,
     $                   MAXIT, NEED1, NEED2
      REAL               TOL, ALPHA, BETA, BNRM2, RHO, RHO1, RHOTOL, 
     $                   GETBREAK, SDOT, SNRM2
*     ..
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     .. External Functions ..
      EXTERNAL           GETBREAK, SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      PHAT = 4
      Q    = 5
      QHAT = 6
      U    = 6
      UHAT = 7
      VHAT = 7
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((QHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((U - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((UHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.9 ) THEN
            NEED1 = ((VHAT - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((QHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((U - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((UHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.9 ) THEN
            NEED2 = ((VHAT - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set breakdown tolerance parameter.
*
      RHOTOL = GETBREAK()
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL SCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R    - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
      ENDIF
*
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
*     Choose RTLD such that initially, (R,RTLD) = RHO is not equal to 0.
*     Here we choose RTLD = R.
*
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform Conjugate Gradient Squared iteration.
*
         ITER = ITER + 1
*
         RHO = SDOT( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
         IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*
*        Compute direction vectors U and P.
*
         IF ( ITER.GT.1 ) THEN
*
*           Compute U.
*
            BETA = RHO / RHO1
            CALL SCOPY( N, WORK(1,R), 1, WORK(1,U), 1 )
            CALL SAXPY( N, BETA, WORK(1,Q), 1, WORK(1,U), 1 )
*
*           Compute P.
*
            CALL SSCAL( N, BETA**2, WORK(1,P), 1 )
            CALL SAXPY( N, BETA, WORK(1,Q), 1, WORK(1,P), 1 )
            CALL SAXPY( N, ONE, WORK(1,U), 1, WORK(1,P), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,R), 1, WORK(1,U), 1 )
            CALL SCOPY( N, WORK(1,U), 1, WORK(1,P), 1 )
         ENDIF
*
*        Compute direction adjusting scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((P    - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,VHAT) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((VHAT - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
         ALPHA = RHO / SDOT( N, WORK(1,RTLD), 1, WORK(1,VHAT), 1 )
*
         CALL SCOPY( N, WORK(1,U), 1, WORK(1,Q), 1 )
         CALL SAXPY( N, -ALPHA, WORK(1,VHAT), 1, WORK(1,Q), 1 )
*
*        Compute direction adjusting vectORT UHAT.
*        PHAT is being used as temporary storage here.
*
         CALL SCOPY( N, WORK(1,Q), 1, WORK(1,PHAT), 1 ) 
         CALL SAXPY( N, ONE, WORK(1,U), 1, WORK(1,PHAT), 1 )
*********CALL PSOLVE( WORK(1,UHAT), WORK(1,PHAT) )
*
         NDX1 = ((UHAT - 1) * LDW) + 1
         NDX2 = ((PHAT - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
*
*****************
 5       CONTINUE
*****************
*
*        Compute new solution approximation vector X.
*
         CALL SAXPY( N, ALPHA, WORK(1,UHAT), 1, X, 1 )
*
*        Compute residual R and check for tolerance.
*
*********CALL MATVEC( ONE, WORK(1,UHAT), ZERO, WORK(1,QHAT) )
*
         NDX1 = ((UHAT - 1) * LDW) + 1
         NDX2 = ((QHAT - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 6
         IJOB = 1
         RETURN
*
*****************
 6       CONTINUE
*****************
*
         CALL SAXPY( N, -ALPHA, WORK(1,QHAT), 1, WORK(1,R), 1 )
*
*********RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 7
         IJOB = 4
         RETURN
*
*****************
 7       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         RHO1 = RHO
*
      GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) INFO = -10
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of CGSREVCOM
*
      END





*
      SUBROUTINE CHEBY( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                  PSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*
*  Purpose
*  =======
*
*  CHEBY solves the linear system Ax = b using the Chebyshev iterative 
*  method with preconditioning. This version requires explicit 
*  knowledge of the maximum and minimum eigenvalues. Note that these
*  eigenvalues must be REAL             positive, which is the case for the
*  symmetric positive definite system.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,3).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*  ============================================================
*     ..
*
*This variable used to communicate requests between CHEBY() and 
*  CHEBYSREVCOM()
*CHEBY -> CHEBYREVCOM: 1 = init, 
*                      2 = use saved state to resume flow.
*CHEBYREVCOM -> CHEBY: -1 = done, return to main, 
*                       1 = matvec using SCLR1/2, NDX1/2 
*                       2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from CHEBYREVCOM to CHEBY
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         CHEBYREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == Z; 3 == P; 4 == Q; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*
*     First time call always init.
      IJOB = 1

 1    CONTINUE

          CALL CHEBYREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call solve.
            CALL PSOLVE(WORK(NDX1), WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 4) THEN
*            do stopping test 2
*           if first time, set INFO so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
            FTFLG = .FALSE.
          ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
          GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of CHEBY.f
*
      END
*
      SUBROUTINE CHEBYREVCOM( N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                        NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( LDW,*)
*     ..
*
*  Purpose
*  =======
*
*  CHEBY solves the linear system Ax = b using the Chebyshev iterative 
*  method with preconditioning. This version requires explicit 
*  knowledge of the maximum and minimum eigenvalues. Note that these
*  eigenvalues must be REAL             positive, which is the case for the
*  symmetric positive definite system.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,3).
*          Workspace for residual, direction vector, etc.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SNRM2
*  ============================================================
*     ..
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, P, Z, MAXIT, NEED1, NEED2
      REAL               ALPHA, BETA, C, D, BNRM2, EIGMAX, EIGMIN, TOL,
     $                   SNRM2
*     ..
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     .. External Routines ..
      EXTERNAL           SAXPY, SCOPY, SNRM2
*
*     .. Executable Statements ..
*
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
*
      MAXIT = ITER
      TOL   = RESID
*
*     Get extremal eigenvalues.
*
      CALL GETEIG( WORK, LDW, EIGMAX, EIGMIN )
*
*     Alias workspace columns.
*
      R = 1
      P = 2
      Z = 3
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using P as temp. storage.
*********CALL SCOPY(N, X, 1, WORK(1,P), 1)
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
         SCLR1 = -ONE
         SCLR2 = ONE
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 30
      ENDIF
*
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
*     Initialize ellipse parameters.
*
      C = ( EIGMAX - EIGMIN ) / TWO
      D = ( EIGMAX + EIGMIN ) / TWO
*  
      ITER = 0
*
   10 CONTINUE
*
*     Perform Chebyshev iteration.
*
         ITER = ITER + 1
*
*********CALL PSOLVE( WORK(1,Z), WORK(1,R) )
*
         NDX1 = ((Z - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
         IF ( ITER.GT.1 ) THEN
            BETA = ( ( C * ALPHA ) / TWO )**2
            ALPHA = ONE / ( D - BETA )
            CALL SAXPY( N, BETA, WORK(1,P), 1, WORK(1,Z), 1 )
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,Z), 1, WORK(1,P), 1 )
            ALPHA = TWO / D
         ENDIF 
*
*        Compute new approximation vector X; check accuracy.
*
         CALL SAXPY( N, ALPHA, WORK(1,P), 1, X, 1 )
*********CALL MATVEC( -ALPHA, WORK(1,P), ONE, WORK(1,R) )
*
         NDX1 = ((P - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = -ALPHA
         SCLR2 = ONE
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
*********RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 5
         IJOB = 4
         RETURN
*
*****************
 5       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of CHEBYREVCOM
*
      END
*
      SUBROUTINE GMRES( N, B, X, RESTRT, WORK, LDW, WORK2, LDW2, 
     $                  ITER, RESID, MATVEC, PSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     EiITERkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, RESTRT, LDW, LDW2, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( * ), WORK2( * )
*     ..
*     .. Function Arguments ..
*
      EXTERNAL           MATVEC, PSOLVE
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix WORK2.
*
*  WORK    (workspace) REAL             array, dimension (LDW,5).
*          Note that if the initial guess is the zero vector, then 
*          storing the initial residual is not necessary.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  WORK2   (workspace) REAL             array, dimension (LDW2,2*RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDW2    (input) INTEGER
*          The leading dimension of the array WORK2.
*          LDW2 >= max(1,RESTRT).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable error tolerance.
*          On ouput, the norm of the residual vector if solution
*          approximated to tolerance, otherwise reset to input
*          tolerance.
*
*  INFO    (output) INTEGER
*          =  0:  successful exit
*          =  1:  maximum number of iterations performed;
*                 convergence not achieved.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product A*x = y.
*          Vector x must remain unchanged. The solution is
*          over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( X, Y )
*
*  ============================================================
*     ..
*     .. Parameters ..
      INTEGER            OFSET
      PARAMETER        ( OFSET = 1000 )
*
*This variable used to communicate requests between GMRES & GMRESREVCOM
*GMRES -> GMRESREVCOM: 1 = init, 
*                      2 = use saved state to resume flow.
*GMRESREVCOM -> GMRES: -1 = done, return to main, 
*                       1 = matvec using SCLR1/2, NDX1/2 
*                       2 = solve using NDX1/2
*                       3 = matvec using WORK2 NDX1 & WORK NDX2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from GMRESREVCOM to GMRES
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         GMRESREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ELSE IF ( LDW2.LT.RESTRT ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*     Note: V, and GIV contain # of history vectors each.
*     To access the i'th vector in V, add i to V*OFSET 1<=i<=RESTRT
*     To access the i'th vector in GIV, add i to GIV*OFSET 1<=i<=RESTRT
*
*     1 == R; 2 == S; 3 == W; 4 == Y; 5 == AV; 6 == H; 7*OFSET+i == V;
*     8*OFSET+i == GIV; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL GMRESREVCOM(N, B, X, RESTRT, WORK, LDW, WORK2, LDW2, 
     $                     ITER, RESID, INFO, NDX1, NDX2, SCLR1, 
     $                     SCLR2, IJOB)


*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec directly with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call solve.
            CALL PSOLVE(WORK(NDX1), WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 4) THEN
*           do stopping test 2
*           if first time, set INFO so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
            FTFLG = .FALSE.
         ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
          GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of GMRES.f
*
      END
*
*
       SUBROUTINE GMRESREVCOM(N, B, X, RESTRT, WORK, LDW, WORK2,
     $                  LDW2, ITER, RESID, INFO, NDX1, NDX2, SCLR1, 
     $                  SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the
*     Solution of Linear Systems: Building Blocks for Iterative
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
*     EiITERkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, RESTRT, LDW, LDW2, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( LDW,* ), WORK2( LDW2,* )
*     ..
*
*  Purpose
*  =======
*
*  GMRES solves the linear system Ax = b using the
*  Generalized Minimal Residual iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*  RESTRT  (input) INTEGER
*          Restart parameter, <= N. This parameter controls the amount
*          of memory required for matrix WORK2.
*
*  WORK    (workspace) REAL             array, dimension (LDW,5).
*          Note that if the initial guess is the zero vector, then 
*          storing the initial residual is not necessary.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  WORK2   (workspace) REAL             array, dimension (LDW2,2*RESTRT+2).
*          This workspace is used for constructing and storing the
*          upper Hessenberg matrix. The two extra columns are used to
*          store the Givens rotation matrices.
*
*  LDW2    (input) INTEGER
*          The leading dimension of the array WORK2.
*          LDW2 >= max(1,RESTRT).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable error tolerance.
*          On ouput, the norm of the residual vector if solution
*          approximated to tolerance, otherwise reset to input
*          tolerance.
*
*  INFO    (output) INTEGER
*          =  0:  successful exit
*          =  1:  maximum number of iterations performed;
*                 convergence not achieved.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  ============================================================
*
*     .. Parameters ..
      REAL                ZERO, ONE
      PARAMETER         ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER             OFSET
      PARAMETER         ( OFSET = 1000 )
*     ..
*     .. Local Scalars ..
      INTEGER             I, MAXIT, AV, GIV, H, R, S, V, W, Y,
     $                    NEED1, NEED2
      REAL                BNRM2, RNORM, TOL, APPROXRES, SDOT, SNRM2
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     ..
*     .. External Routines ..
      EXTERNAL            SAXPY, SCOPY, SDOT, SNRM2
*     ..
*     .. Executable Statements ..
*
* Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 200
      ENDIF
*
* init.
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R   = 1
      S   = 2
      W   = 3
      Y   = 4
      AV  = 5
      V   = 6
*
      H   = 1
      GIV = H + RESTRT
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((W - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((AV - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( ( NDX1.GT.V*OFSET ) .AND.
     $           ( NDX1.LE.V*OFSET+RESTRT ) ) THEN
            NEED1 = ((NDX1-V*OFSET - 1) * LDW) + 1
         ELSEIF( ( NDX1.GT.GIV*OFSET ) .AND.
     $           ( NDX1.LE.GIV*OFSET+RESTRT ) ) THEN
            NEED1 = ((NDX1-GIV*OFSET - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 100
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((W - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((AV - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( ( NDX2.GT.V*OFSET ) .AND.
     $           ( NDX2.LE.V*OFSET+RESTRT ) ) THEN
            NEED2 = ((NDX2-V*OFSET - 1) * LDW) + 1
         ELSEIF( ( NDX2.GT.GIV*OFSET ) .AND.
     $           ( NDX2.LE.GIV*OFSET+RESTRT ) ) THEN
            NEED2 = ((NDX2-GIV*OFSET - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 100
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using X directly
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 1
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ).LT.TOL ) GO TO 200
      ENDIF
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
      ITER = 0
   10 CONTINUE
*
         ITER = ITER + 1
*
*        Construct the first column of V, and initialize S to the
*        elementary vector E1 scaled by RNORM.
*
*********CALL PSOLVE( WORK( 1,V ), WORK( 1,R ) )
*
         NDX1 = ((V - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         RNORM = SNRM2( N, WORK( 1,V ), 1 )
         CALL SSCAL( N, RNORM, WORK( 1,V ), 1 )
         CALL ELEMVEC( 1, N, RNORM, WORK( 1,S ) )
*
         DO 50 I = 1, RESTRT
************CALL MATVEC( ONE, WORK( 1,V+I-1 ), ZERO, WORK( 1,AV ) )
*
         NDX1 = ((V+I-1 - 1) * LDW) + 1
         NDX2 = ((AV    - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 3
         RETURN
*
*****************
 4       CONTINUE
*****************
*
*********CALL PSOLVE( WORK( 1,W ), WORK( 1,AV ) )
*
         NDX1 = ((W  - 1) * LDW) + 1
         NDX2 = ((AV - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 5
         IJOB = 2
         RETURN
*
*****************
 5       CONTINUE
*****************
*
*           Construct I-th column of H so that it is orthnormal to 
*           the previous I-1 columns.
*
            CALL ORTHOH( I, N, WORK2( 1,I+H-1 ), WORK( 1,V ), LDW,
     $                   WORK( 1,W ) )
*
            IF ( I.GT.1 )
*
*              Apply Givens rotations to the I-th column of H. This
*              effectively reduces the Hessenberg matrix to upper
*              triangular form during the RESTRT iterations.
*
     $         CALL APPLYGIVENS( I, WORK2( 1,I+H-1 ), WORK2( 1,GIV ),
     $                           LDW2 )
*
*           Approxiyate residual norm. Check tolerance. If okay, compute
*           final approximation vector X and quit.
*
            RESID = APPROXRES( I, WORK2( 1,I+H-1 ), WORK( 1,S ),
     $                         WORK2( 1,GIV ), LDW2 ) / BNRM2
            IF ( RESID.LE.TOL ) THEN
               CALL UPDATE(I, N, X, WORK2( 1,H ), LDW2, 
     $                     WORK(1,Y), WORK(1,S), WORK( 1,V ), LDW)
               GO TO 200
            ENDIF
   50    CONTINUE
*
*        Compute current solution vector X.
*
         CALL UPDATE(RESTRT, N, X, WORK2( 1,H ), LDW2,
     $               WORK(1,Y), WORK( 1,S ), WORK( 1,V ), LDW )
*
*        Compute residual vector R, find norm,
*        then check for tolerance.
*
         CALL SCOPY( N, B, 1, WORK( 1,R ), 1 )
*********CALL MATVEC( -ONE, X, ONE, WORK( 1,R ) )
*
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = -ONE
         SCLR2 = ONE
         RLBL = 6
         IJOB = 1
         RETURN
*
*****************
 6       CONTINUE
*****************
*
         WORK( I+1,S ) = SNRM2( N, WORK( 1,R ), 1 )
*
*********RESID = WORK( I+1,S ) / BNRM2
*********IF ( RESID.LE.TOL  ) GO TO 200
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 7
         IJOB = 4
         RETURN
*
*****************
 7       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 200
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 100
         ENDIF
*
         GO TO 10
*
  100 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
  200 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1

      RETURN
*
*     End of GMRESREVCOM
*
      END
*
*     =========================================================
      SUBROUTINE ORTHOH( I, N, H, V, LDV, W )
*
      INTEGER            I, N, LDV
      REAL               H( * ), W( * ), V( LDV,* )
*
*     Construct the I-th column of the upper Hessenberg matrix H
*     using the Gram-Schmidt process on V and W.
*
      INTEGER            K
      REAL               SDOT, SNRM2, ONE
      PARAMETER        ( ONE = 1.0E+0 )
      EXTERNAL           SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*
      DO 10 K = 1, I
         H( K ) = SDOT( N, W, 1, V( 1,K ), 1 )
         CALL SAXPY( N, -H( K ), V( 1,K ), 1, W, 1 )
   10 CONTINUE
      H( I+1 ) = SNRM2( N, W, 1 )
      CALL SCOPY( N, W, 1, V( 1,I+1 ), 1 )
      CALL SSCAL( N, ONE / H( I+1 ), V( 1,I+1 ), 1 )
*
      RETURN
*
      END
*     =========================================================
      SUBROUTINE APPLYGIVENS( I, H, GIVENS, LDG )
*
      INTEGER            I, LDG
      REAL               H( * ), GIVENS( LDG,* )
*
*     This routine applies a sequence of I-1 Givens rotations to
*     the I-th column of H. The Givens parameters are stored, so that
*     the first I-2 Givens rotation matrices are known. The I-1st
*     Givens rotation is computed using BLAS 1 routine SROTG. Each
*     rotation is applied to the 2x1 vector [H( J ), H( J+1 )]',
*     which results in H( J+1 ) = 0.
*
      INTEGER            J
*      REAL               TEMP
      EXTERNAL           SROTG
*
*     .. Executable Statements ..
*
*     Construct I-1st rotation matrix.
*
*     CALL SROTG( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL GETGIV( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*
*     Apply 1,...,I-1st rotation matrices to the I-th column of H.
*
      DO 10 J = 1, I-1
         CALL ROTVEC( H( J ), H( J+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*        TEMP     =  GIVENS( J,1 ) * H( J ) + GIVENS( J,2 ) * H( J+1 ) 
*        H( J+1 ) = -GIVENS( J,2 ) * H( J ) + GIVENS( J,1 ) * H( J+1 )
*        H( J ) = TEMP
 10   CONTINUE
*
      RETURN
*
      END
*
*     ===============================================================
      REAL             FUNCTION APPROXRES( I, H, S, GIVENS, LDG )
*
      INTEGER            I, LDG
      REAL               H( * ), S( * ), GIVENS( LDG,* )
*
*     This function allows the user to approximate the residual
*     using an updating scheme involving Givens rotations. The
*     rotation matrix is formed using [H( I ),H( I+1 )]' with the
*     intent of zeroing H( I+1 ), but here is applied to the 2x1
*     vector [S(I), S(I+1)]'.
*
      INTRINSIC          ABS
      EXTERNAL           SROTG
*
*     .. Executable Statements ..
*
*     CALL SROTG( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL GETGIV( H( I ), H( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
      CALL ROTVEC( S( I ), S( I+1 ), GIVENS( I,1 ), GIVENS( I,2 ) )
*
      APPROXRES = ABS( S( I+1 ) )
*
      RETURN
*
      END
*     ===============================================================
      SUBROUTINE UPDATE( I, N, X, H, LDH, Y, S, V, LDV )
*
      INTEGER            N, I, J, LDH, LDV
      REAL               X( * ), Y( * ), S( * ), H( LDH,* ), V( LDV,* )
      EXTERNAL           SAXPY, SCOPY, STRSV
*
*     Solve H*y = s for upper triangualar H.
*
      CALL SCOPY( I, S, 1, Y, 1 )
      CALL STRSV( 'UPPER', 'NOTRANS', 'NONUNIT', I, H, LDH, Y, 1 )
*
*     Compute current solution vector X.
*
      DO 10 J = 1, I
         CALL SAXPY( N, Y( J ), V( 1,J ), 1, X, 1 )
   10 CONTINUE
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE GETGIV( A, B, C, S )
*
      REAL               A, B, C, S, TEMP, ZERO, ONE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*
      IF ( B.EQ.ZERO ) THEN
         C = ONE
         S = ZERO
      ELSE IF ( ABS( B ).GT.ABS( A ) ) THEN
         TEMP = -A / B
         S = ONE / SQRT( ONE + TEMP**2 )
         C = TEMP * S
      ELSE
         TEMP = -B / A
         C = ONE / SQRT( ONE + TEMP**2 )
         S = TEMP * C
      ENDIF
*
      RETURN
*
      END
*
*     ================================================================
      SUBROUTINE ROTVEC( X, Y, C, S )
*
      REAL               X, Y, C, S, TEMP

*
      TEMP = C * X - S * Y
      Y    = S * X + C * Y
      X    = TEMP
*
      RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE ELEMVEC( I, N, ALPHA, E )
*
*     Construct the I-th elementary vector E, scaled by ALPHA.
*
      INTEGER            I, J, N
      REAL               ALPHA, E( * )
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*
      DO 10 J = 1, N
         E( J ) = ZERO
   10 CONTINUE
      E( I ) = ALPHA
*
      RETURN
*
      END
*
      SUBROUTINE JACOBI( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                   INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC
*
*  Purpose
*  =======
*
*  JACOBI solves the linear system Ax = b using the Jacobi iterative 
*  method. The matrix splitting should be accomplished before calling
*  this routine. The diagonal elements of the matrix must be passed into
*  this routine in the first column of matrix WORK. 
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,3)
*          Workspace for residual, direction vector, etc.
*          
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( x - x_1 ) / norm( x ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*  ==========================================================
*
*     .. Local Scalars ..
*
*This variable used to communicate requests between JACOBI() &
* JACOBIREVCOM()
*JACOBI -> JACOBIREVCOM: 1 = init, 
*                        2 = use saved state to resume flow.
*JACOBIREVCOM -> JACOBI: -1 = done, return to main, 
*                         1 = matvec using SCLR1/2, NDX1/2 
      INTEGER          IJOB
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from JacobiREVCOM to Jacobi.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, XNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         JACOBIREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == MM; 2 == X1; 3 == TEMP; -1 == ignore; any other == error
      NDX1 = 2
      NDX2 = -1
      TOL = RESID
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL JACOBIREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                      NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
*           note: ndx1 == -1 here, as we're directly using X
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 2) THEN
*           do stopping test 2
*           XNRM2 shld be recomputed everytime.
            INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), X, XNRM2, RESID, TOL, INFO)
         ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
      RETURN
*
*     End of JACOBI
*
      END
*
      SUBROUTINE JACOBIREVCOM(N, B, X, WORK, LDW, ITER, RESID,
     $                     INFO, NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( LDW,* )
*
*  Purpose
*  =======
*
*  JACOBI solves the linear system Ax = b using the Jacobi iterative 
*  method. The matrix splitting should be accomplished before calling
*  this routine. The diagonal elements of the matrix must be passed into
*  this routine in the first column of matrix WORK. 
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) REAL             array, dimension (LDW,3)
*          Workspace for residual, direction vector, etc.
*          
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( x - x_1 ) / norm( x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SNRM2
*  ==========================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, MAXIT, X1, MM, TEMP, NEED1, NEED2
      REAL               TOL, SNRM2
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Routines ..
*
      EXTERNAL           SAXPY, SCOPY, SNRM2, MATSPLIT
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      MM   = 1
      X1   = 2
      TEMP = 3
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((MM - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((TEMP - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((MM - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((TEMP - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
      ITER = 0
*
*     Form matrix splitting inv(M) and N.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW,'JACOBI','SPLIT' )
*
   10 CONTINUE
*
*        Perform Jacobi iteration
*
         ITER = ITER + 1
*
*        Save the current approximation to X in X1.
*
         CALL SCOPY( N, X, 1, WORK( 1,X1 ),  1 )
*
*        Apply iteration; result is updated approximation vector x.
*
         CALL SCOPY( N, B, 1, WORK( 1,TEMP ), 1 )
*********CALL MATVEC( ONE, X, ONE, WORK( 1,TEMP ) )
*
         NDX1 = -1
         NDX2 = ((TEMP - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ONE
         RLBL = 2
         IJOB = 1
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         DO 15 I = 1, N
            X( I ) = WORK( I,MM ) *  WORK( I,TEMP )
   15    CONTINUE
*
*        Compute error and check for acceptable convergence.
*
         CALL SAXPY( N, -ONE, X, 1, WORK( 1,X1 ), 1 )
*
*********RESID = SNRM2( N, WORK( 1,X1 ), 1 ) / SNRM2( N, X, 1 )
*********IF ( RESID.LE.TOL  ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
      GO TO 10
*
   20 CONTINUE
*
*     Iteration fails.
*     Reconstruct matrix A.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW, 'JACOBI','RECONSTRUCT' )
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   30 CONTINUE 
*
*     Iteration successful. Reconstruct matrix A.
*
      CALL MATSPLIT( ONE, B, WORK( 1,MM ), LDW, 'JACOBI','RECONSTRUCT' )
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of JACOBIREVCOM
*
      END


*
*     This file contains routines for forming matrices that result from
*     a 5- or 7-point discretization of elliptic PDEs with Dirichlet
*     boundary conditions, and a consistent mass matrix "Wathen".
*
*        GEN57PT and GETSTEN are from SPARSEKIT. They actually form the
*        row compressed matrix.
*
*        COEFF provides the functions for computing the coefficients 
*        of the PDE.
*
*        Finally, for testing the iterative templates, COMP2DENSE converts
*        the row compressed matrix to dense form.
*
*     =================================================================
      SUBROUTINE GEN57PT( NX, NY, NZ, A, INDX, POINTR, AFUN, BFUN, CFUN,
     $                    DFUN, EFUN, FFUN, GFUN )
*
*     .. Scalar Arguments ..
      INTEGER            NX, NY, NZ
*     ..
*     .. Array Arguments ..
      INTEGER            INDX( * ), POINTR( * )
      REAL               A( * )
*     ..
*     .. External Functions ..
      EXTERNAL           AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN
*
*  Purpose
*  =======
*
*  Adapted/altered from SPARSEKIT
*
*  This subroutine computes the sparse matrix in row compressed
*  format for the elliptic operator
*
*  L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
*        delx ( d u ) + dely (e u) + delz( f u ) + g u 
*
*  with Dirichlet Boundary conditions, on a rectangular 1-D, 
*  2-D or 3-D grid using centered difference schemes.
* 
*  The functions a, b, ..., g are known through the
*  subroutines  afun, bfun, ..., gfun.
*  Note that to obtain the correct matrix, any function that is not
*  needed should be set to zero. For example for two-dimensional
*  problems, nz should be set to 1 and the functions cfun and ffun
*  should be zero functions. 
*
*  Uses natural ordering, first x direction, then y, then z
*  mesh size h is uniform and determined by grid points 
*  in the x-direction.
*
*  Arguments
*  =========
*
*  NX     (input) INTEGER
*         Number of points in X direction.
*
*  NY     (input) INTEGER
*         Number of points in Y direction.
*
*  NZ     (input) INTEGER
*         Number of points in Z direction.
*
*  A,     (output) REAL array.
*         Nonzero elements of the matrix. Stored in row compressed form.
* 
*  INDX   (output) INTEGER array.
*         Column index of matrix element.
* 
*  POINTR (output) INTEGER array.
*         Each element = P+1, where P is the number of nonzero elements
*         in the preceding rows of the matrix.
*
*  AFUN,
*  BFUN,
*  CFUN,
*  DFUN,
*  EFUN,
*  FFUN,
*  GFUN   (external subroutine) 
*         The user must supply the functions for computing the coefficients
*         of the PDE.
*
*  Description of the STENCIL:
*
*     stencil [1:7] has the following meaning:
*
*        center point = stencil(1)
*        west point   = stencil(2)
*        east point   = stencil(3)
*        south point  = stencil(4)
*        north point  = stencil(5)
*        front point  = stencil(6) 
*        back point   = stencil(7)
*
*                           st(5)
*                            |
*                            |  
*                            |
*                            |          .st(7)
*                            |     .
*                            | . 
*         st(2) ----------- st(1) ---------- st(3)
*                       .    |
*                   .        |
*               .            |
*            st(6)           |
*                            |
*                            |
*                           st(4)
*
*     ===============================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            KX, KY, KZ, IX, IY, IZ, IEDGE, NODE
      REAL               H, STENCIL( 7 )
*
*     .. Executable Statements ..
*
*     Initializations
*
      H     = ONE / (NX+1)
      KX    = 1
      KY    = NX
      KZ    = NX*NY
      IEDGE = 1
      NODE  = 1
*
      DO 30 IZ = 1, NZ
         DO 20 IY = 1, NY
            DO 10 IX = 1, NX
               POINTR(NODE) = IEDGE
*
*              Get stencil.
*
               CALL GETSTEN( NY, NZ, IX, IY, IZ, STENCIL, H,
     $                       AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN )
*
*              West.
*
               IF ( IX.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KX
                  A(IEDGE) = STENCIL(2)
                  IEDGE=IEDGE + 1
               ENDIF
*
*              South.
*
               IF ( IY.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KY
                  A(IEDGE) = STENCIL(4)
                  IEDGE=IEDGE + 1
               END IF
*
*              Front Plane.
*
               IF ( IZ.GT.1 ) THEN
                  INDX(IEDGE)=NODE-KZ
                  A(IEDGE) = STENCIL(6)
                  IEDGE=IEDGE + 1
               ENDIF
*
*              Center node.
*
               INDX(IEDGE) = NODE
               A(IEDGE) = STENCIL(1)
               IEDGE = IEDGE + 1
*
*              Upper part.
*
*              East.
*
               IF ( IX.LT.NX ) THEN
                  INDX(IEDGE)=NODE+KX
                  A(IEDGE) = STENCIL(3)
                  IEDGE=IEDGE + 1
               END IF
               IF ( IY.LT.NY ) THEN
                  INDX(IEDGE)=NODE+KY
                  A(IEDGE) = STENCIL(5)
                  IEDGE=IEDGE + 1
               END IF
*
*              Back plane.
*
               IF ( IZ.LT.NZ ) THEN
                  INDX(IEDGE)=NODE+KZ
                  A(IEDGE) = STENCIL(7)
                  IEDGE=IEDGE + 1
               END IF
*
*              Next node.
*
               NODE=NODE+1
*
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
*
      POINTR(NODE)=IEDGE
*
      RETURN
*
*     -- End of GEN57PT
*
      END
*     ===============================================================
      SUBROUTINE GETSTEN( NY, NZ, KX, KY, KZ, STENCIL, H,
     $                    AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN )
*
*     .. Argument Declarations ..
      INTEGER            NY, NZ, KX, KY, KZ
      REAL               H, AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN,
     $                   STENCIL( * )
*     ..
*     .. External Functions ..
      EXTERNAL           AFUN, BFUN, CFUN, DFUN, EFUN, FFUN, GFUN
*
*  Purpose
*  =======
*
*  This subroutine calcultes the correct stencil values for
*  elliptic operator
*
*     L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) + 
*           delx ( d u ) + dely (e u) + delz( f u ) + g u.
*
*  For 2-D problems the discretization formula that is used is:
*      
*  h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
*               a(i-1/2,j)*{u(i-1,j) - u(i,j)} + 
*               b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
*               b(i,j-1/2)*{u(i,j-1) - u(i,j)} + 
*              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
*              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} + 
*              (h**2)*g(i,j)*u(i,j) 
*
*  ===================================================================
*
*     .. Parameters ..
      REAL               ZERO, HALF
      PARAMETER        ( ZERO = 0.0E+0, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            K
      REAL               HHALF,CNTR, X, Y, Z, COEFF
*
*     .. Executable Statements ..
*
      DO 10 K=1,7
         STENCIL(K) = ZERO
  10  CONTINUE
*
      HHALF = H*HALF
      X = H*KX
      Y = H*KY
      Z = H*KZ
      CNTR = ZERO
*
*     Differentiation w.r.t. X.
*
      COEFF = AFUN( X+HHALF, Y, Z )
      STENCIL(3) = STENCIL(3) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = AFUN( X-HHALF,Y,Z )
      STENCIL(2) = STENCIL(2) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = DFUN( X, Y, Z )*HHALF
      STENCIL(3) = STENCIL(3) + COEFF
      STENCIL(2) = STENCIL(2) - COEFF
      IF (NY.LE.1) GOTO 99
*
*     Differentiation w.r.t. Y.
*
      COEFF = BFUN( X, Y+HHALF, Z )
      STENCIL(5) = STENCIL(5) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = BFUN( X, Y-HHALF, Z )
      STENCIL(4) = STENCIL(4) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = EFUN( X, Y, Z )*HHALF
      STENCIL(5) = STENCIL(5) + COEFF
      STENCIL(4) = STENCIL(4) - COEFF
      IF ( NZ.LE.1) GOTO 99
*
*     Differentiation w.r.t. Z.
*
      COEFF = CFUN( X, Y, Z+HHALF )
      STENCIL(7) = STENCIL(7) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = CFUN( X, Y, Z-HHALF )
      STENCIL(6) = STENCIL(6) + COEFF
      CNTR = CNTR + COEFF
*
      COEFF = FFUN( X, Y,Z )*HHALF
      STENCIL(7) = STENCIL(7) + COEFF
      STENCIL(6) = STENCIL(6) - COEFF
*
*     Discretization of product by G.
*
 99   COEFF = GFUN( X, Y, Z )
      STENCIL(1) = H*H*COEFF - CNTR
*
      RETURN
*
      END
*     =============================================================
*     Below are some functions for computing the value of the
*     coefficients.
*     =============================================================
      REAL             FUNCTION ZEROFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Function to return ZERO.
*     =======
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*
*     .. Executable Statements ..
*
      ZEROFUN = ZERO
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION ONEFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Function to return ONE.
*     =======
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE = 1.0E+0 )
*
*     .. Executable Statements ..
*
      ONEFUN = ONE
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION NEGONEFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Function to return -ONE.
*     =======
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE = 1.0E+0 )
*
*     .. Executable Statements ..
*
      NEGONEFUN = -ONE
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION THOUSFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Function to return 1000.
*     =======
*
*     .. Executable Statements ..
*
      THOUSFUN = 1000.0E+0
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION TEN5X2FUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Function to return 10 * X^2.
*     =======
*
*     .. Executable Statements ..
*
      TEN5X2FUN = 10.0E+5 * X * X
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION THOUSXFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
* 
*     Purpose: Evaluates the coefficient function
*     =======
*
*     .. Parameter ..
      REAL               NATX
      PARAMETER        ( NATX = 2.71828174E+0 )
*
*     .. Executable Statements ..
*
      THOUSXFUN = 1000.E+0 * ( NATX**( X*Y*Z ) )
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION NEGTHOUSXFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Evaluates the coefficient function
*     =======
*
*     .. Parameter ..
      REAL               NATX
      PARAMETER        ( NATX = 2.71828174E+0 )
*
*     .. Executable Statements ..
*
      NEGTHOUSXFUN = -1000.E+0 * ( NATX**( X*Y*Z ) )
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION HENKFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Evaluates the derivative of the above coefficient function
*     =======
*
*     .. Parameter ..
      REAL               NATX
      PARAMETER        ( NATX = 2.71828174E+0 )
*
*     .. Executable Statements ..
*
      HENKFUN = 20 * ( NATX**( 3.5E+0 * ( X**2 + Y**2 ) ) )
*
*     RETURN
*
      END
*
*     =============================================================
      REAL             FUNCTION HENKDFUN( X, Y, Z )
*
*     .. Argument Declarations ..
      REAL               X, Y, Z
*
*     Purpose: Evaluates the derivative of the above coefficient function.
*     =======
*
*     .. Parameter ..
      REAL               NATX
      PARAMETER        ( NATX = 2.71828174E+0 )
*
*     .. Executable Statements ..
*
      HENKDFUN = 70 * X * ( NATX**( 3.5E+0 * ( X**2 + Y**2 ) ) )
*
*     RETURN
*
      END
*
*     ===============================================================
      SUBROUTINE COMP2DENSE( ASPARSE, POINTR, INDX, N, ADENSE, LDA,
     $                       FLAG, INFO  )
*
      INTEGER            N, LDA, INFO, POINTR( * ), INDX( * )
      REAL               ASPARSE( * ), ADENSE( LDA,* )
      CHARACTER          FLAG*3
*
*     Convert sparse matrix storage to dense.
*
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
      INTEGER            I, J
      LOGICAL            LSAME, LSAMEN
*
*     .. Executable Statements ..
*
      INFO = 0
      IF ( N.LE.0 ) THEN
         INFO = -1
      ELSE IF ( N.GT.LDA ) THEN
         INFO = -2
      ELSE IF ( .NOT.LSAMEN( 3, FLAG,'ROW' ).OR.
     $          .NOT.LSAMEN( 3, FLAG,'ROW' ) ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            ADENSE( I,J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
      IF ( LSAME( FLAG,'ROW' ) ) THEN
         DO 40 I = 1, N
            DO 30 J = POINTR( I ), POINTR( I+1 ) - 1
               ADENSE( I,INDX( J ) ) = ASPARSE( J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF ( LSAME( FLAG,'COL' ) ) THEN
         DO 60 J = 1, N
            DO 50 I = POINTR( J ), POINTR( J+1 ) - 1
               ADENSE( INDX( I ),J ) = ASPARSE( I ) 
   50       CONTINUE
   60    CONTINUE
      ENDIF
*
      RETURN
*
      END
*
*     ================================================================
      SUBROUTINE WATHEN( NX, NY, KK, N, A, LDA, WORK, LDW, INFO )
*
      INTEGER            NX, NY, KK, N, LDA, LDW, INFO
      REAL               A( LDA,* ), WORK( LDW,* )
*
*     Translated from the matlab version found on netlib.
*
*     A is a random N-by-N finite element matrix where 
*     N = 3*NX*NY + 2*NX + 2*NY + 1. A is precisely the "consistent 
*     mass matrix" for a regular NX-by-NY grid of 8-node (serendipity) 
*     elements in 2 space dimensions. A is symmetric positive definite 
*     for any (positive) values of the "density", RHO(NX,NY), which is 
*     chosen randomly in this routine. In particular, if D=DIAG(DIAG(A)), 
*     then 0.25 <= EIG(INV(D)*A) <= 4.5 for any positive integers NX and NY 
*     and any densities RHO(NX,NY). This diagonally scaled matrix is 
*     returned by WATHEN(NX,NY,1).
*
*     Reference: A.J.Wathen, REAListic eigenvalue bounds for 
*     the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987), pp. 449-457.
*
*     BEWARE - this is a sparse matrix, stored in -dense- form, and 
*              it quickly gets large!
*
*     .. Local Scalars ..
*
      INTEGER            I, J, E, EM, RHO, KROW, KCOL,
     $                   ISEED( 4 ), NN( 8 )
      REAL               RHOIT, SLARAN, ZERO, ONE
      PARAMETER         ( ZERO = 0.0E+0, ONE = 1.0E+0  )
*
*     .. Executable Statements ..
*
      INFO = 0
      N = 3*NX*NY + 2*NX + 2*NY + 1
      IF ( N.GT.LDA ) THEN
         WRITE(*,*) 'NOT ENOUGH ROOM ALLOCATED FOR WATHEN MATRIX'
         INFO = -1
         RETURN
      ELSE IF ( NX.LT.1 ) THEN
         INFO = -2
      ELSE IF ( NY.LT.1 ) THEN
         INFO = -3
      ELSE IF ( MAX( NX,NY ).GT.LDW ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Alias workspace columns.
*
      E   = 1
      EM  = E + 8
      RHO = EM + 8
*
      CALL SET_E( WORK( 1,E ), LDW )
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            A( I,J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
      ISEED( 1 ) = 304
      ISEED( 2 ) = 152
      ISEED( 3 ) = 2042
      ISEED( 4 ) = 77
      DO 40 J = 1, NY
         DO 30 I = 1, NX
            WORK( I,RHO+J-1 ) = 100*SLARAN( ISEED )
   30    CONTINUE
   40 CONTINUE
*
      DO 100 J = 1, NY
         DO 90 I = 1, NX
*
            NN(1) = 3*J*NX + 2*I + 2*J + 1
            NN(2) = NN(1)-1
            NN(3) = NN(2)-1
            NN(4) = (3*J-1)*NX+2*J+I-1
            NN(5) = 3*(J-1)*NX+2*I+2*J-3
            NN(6) = NN(5)+1
            NN(7) = NN(6)+1
            NN(8) = NN(4)+1
*
            RHOIT = WORK( I,RHO+J-1 )
            DO 60 KROW = 1, 8
               DO 50 KCOL = 1, 8
                  WORK( KROW,EM+KCOL-1 ) = RHOIT * WORK( KROW,E+KCOL-1 )
   50          CONTINUE
   60       CONTINUE
*
            DO 80 KROW = 1, 8
               DO 70 KCOL = 1, 8
                  A(NN(KROW),NN(KCOL)) = A(NN(KROW),NN(KCOL)) +
     $                                      WORK( KROW,EM+KCOL-1)
   70          CONTINUE
   80       CONTINUE 
*
   90    CONTINUE
  100 CONTINUE
*  
      IF ( KK.EQ.1 ) THEN
*
*        A = diag(diag(A)) \ A (the result being unit diagonal);
*
         DO 110 I = 1, J
            A( I,I ) = ONE
  110    CONTINUE 
      ENDIF
*
      RETURN
*
      END
*     ========================================================
      SUBROUTINE SET_E( E, LDE )
*
      INTEGER            I, J, LDE
      REAL               SCALE, E( LDE,* )
*
      EXTERNAL           SSCAL
*
      E( 1,1 ) =  6.0E+0
      E( 2,1 ) = -6.0E+0
      E( 3,1 ) =  2.0E+0
      E( 4,1 ) = -8.0E+0
      E( 1,2 ) = -6.0E+0
      E( 2,2 ) = 32.0E+0
      E( 3,2 ) = -6.0E+0
      E( 4,2 ) = 20.0E+0
      E( 1,3 ) =  2.0E+0
      E( 2,3 ) = -6.0E+0
      E( 3,3 ) =  6.0E+0
      E( 4,3 ) = -6.0E+0
      E( 1,4 ) = -8.0E+0
      E( 2,4 ) = 20.0E+0
      E( 3,4 ) = -6.0E+0
      E( 4,4 ) = 32.0E+0
*
      E( 1,5 ) =  3.0E+0
      E( 2,5 ) = -8.0E+0
      E( 3,5 ) =  2.0E+0
      E( 4,5 ) = -6.0E+0
      E( 1,6 ) = -8.0E+0
      E( 2,6 ) = 16.0E+0
      E( 3,6 ) = -8.0E+0
      E( 4,6 ) = 20.0E+0
      E( 1,7 ) =  2.0E+0
      E( 2,7 ) = -8.0E+0
      E( 3,7 ) =  3.0E+0
      E( 4,7 ) = -8.0E+0
      E( 1,8 ) = -6.0E+0
      E( 2,8 ) = 20.0E+0
      E( 3,8 ) = -8.0E+0
      E( 4,8 ) = 16.0E+0
*
      DO 20 J = 1, 4
         DO 10 I = 5, 8 
            E( I,J ) = E( J,I ) 
   10    CONTINUE
   20 CONTINUE
*
      DO 40 J = 5, 8
         DO 30 I = 5, 8 
            E( I,J ) = E( I-4,J-4 )
   30    CONTINUE
   40 CONTINUE
*
      SCALE = 1.0E+0 / 45.0E+0
      DO 50 I = 1, 8
         CALL SSCAL( 8, SCALE, E( 1,I ), 1 )
   50 CONTINUE
*
      RETURN
*
      END







*
      SUBROUTINE MATVEC( ALPHA, X, BETA, Y )
*
*     This MatVec routine assumes the matrix is in dense format,
*     and uses the BLAS SGEMV.
*
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               X( * ), Y( * )
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      REAL                A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
*     ..
      EXTERNAL           SGEMV
*
      CALL SGEMV( 'NOTRANSPOSE', N, N, ALPHA, A, LDA, X, 1, BETA, Y, 1 )
*
      RETURN
*
      END
*
*     =================================================
      SUBROUTINE MATVECTRANS( ALPHA, X, BETA, Y )
*
*     This MatVec routine assumes the matrix is in dense format,
*     and uses the BLAS SGEMV.
*
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               X( * ), Y( * )
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      REAL                A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
*     ..
      EXTERNAL           SGEMV
*
      CALL SGEMV( 'TRANSPOSE', N, N, ALPHA, A, LDA, X, 1, BETA, Y, 1 )
*
      RETURN
*
      END
*
*  This file contains the preconditioner solve routines:
*
*  PSOLVE and PSOLVETRANS call the appropriate solver:
*
*     PSOLVENONE and PSOLVENONETRANS for using no preconditioning.
*
*     PSOLVEJAC and PSOLVEJACTRANS for Jacobi preconditioning.
*
*     Also included are the solvers for QMR which require left and right
*     preconditioning: PSOLVEQ and PSOLVETRANSQ
*
*     ================================================================
      SUBROUTINE PSOLVE( X, B )
*
*     .. Array Arguments ..
      REAL                X( * ), B( * )
*     ..
*     .. Local Scalars ..
      LOGICAL             LSAME
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         CALL PSOLVEJAC( X, B )
      ELSE
         WRITE(*,*) 'IN PSOLVE: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVETRANS( X, B )
*
*     .. Array Arguments ..
      REAL                X( * ), B( * )
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*     ..
*     .. Local Scalars ..
      LOGICAL          LSAME
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONETRANS( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         CALL PSOLVEJACTRANS( X, B )
      ELSE
         WRITE(*,*) 'IN PSOLVE: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVENONE( X, B )
*
*     .. Array Arguments ..
      REAL             X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             N, LDA
      COMMON            / MATDIM / N, LDA
*
*  Purpose
*  =======
*
*  This PSOLVE is for the unpreconditioned version, i.e. just does
*  a vector copy ( B to X ) then returns.
*
*  Arguments
*  =========
*
*  B       (input) REAL array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) REAL array, dimension N.
*          Set to solution on output.
*
*  BLAS:  SCOPY
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL SCOPY( N, B, 1, X, 1 )
*
      RETURN
*
*     End of PSolveNone
*
      END
*
*     =====================================================
      SUBROUTINE PSOLVENONETRANS( X, B )
*
*     .. Array Arguments ..
*     ..
      REAL             X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             N, LDA
      COMMON            / MATDIM / N, LDA
*     ..
*
*  Purpose
*  =======
*
*  This PSOLVE is for the unpreconditioned version, i.e. just does
*  a vector copy ( B to X ) then returns.
*
*  Arguments
*  =========
*
*  B       (input) REAL array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) REAL array, dimension N.
*          Set to solution on output.
*
*  BLAS:  SCOPY
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL SCOPY( N, B, 1, X, 1 )
*
      RETURN
*
*     End of PSolve
*
      END
*
*     ===========================================================
      SUBROUTINE PSOLVEJAC( X, B )
*
*     .. Array Delcarations ..
      REAL                X( * ), B( * )
*     ..
*     .. Common Blocks ..
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      REAL                A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
*
*  Purpose
*  =======
*
*  PSOLVE solves the linear system Mx = b where matrix M has 
*  is diagonal. Note that we assume M actually equals inv(M),
*  so that we may multiply, rather than divide.
*
*  Arguments
*  =========
*
*  B       (input) REAL array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) REAL array, dimension N.
*          Set to solution on output.
*  ============================================================
*
*     .. Local Scalars ..
*
      INTEGER          I
*
*     .. Executable Statements ..
*
      DO 10 I = 1, N
         X( I ) = B( I ) / M( I )
   10 CONTINUE
*
      RETURN
*
*     End of PSolveJac
*
      END
*
*     =========================================================
      SUBROUTINE PSOLVEJACTRANS( X, B )
*
*     .. Array Delcarations ..
*
      REAL             X( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  PSOLVETRANS solves the linear system Mx = b where matrix M has
*  is diagonal. Since this is the same as the non-transpose version,
*  this routine is actual just a mask to PSOLVEJAC.
*
*  Arguments
*  =========
*
*  B       (input) REAL array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (output) REAL array, dimension N.
*          Set to solution on output.
*  ============================================================
*
*     .. Executable Statements ..
*
      CALL PSOLVEJAC( X, B )
*
      RETURN
*
*     End of PSolveJacTrans
*
      END
*
*     ================================================================
      SUBROUTINE PSOLVEQ( X, B, WHICH )
*
*     .. Argument Declarations ..
      REAL                X( * ), B( * )
      CHARACTER           WHICH*4
*     ..
*     .. Local Scalars ..
      LOGICAL             LSAME
*     ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         IF( LSAME( WHICH, 'LEFT' ) ) THEN
            CALL PSOLVEJAC( X, B )
         ELSE
            CALL PSOLVENONE( X, B )
         ENDIF
      ELSE
         WRITE(*,*) 'IN PSOLVEQ: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*     ================================================================
*
      SUBROUTINE PSOLVETRANSQ( X, B, WHICH )
*
*     .. Argument Declarations ..
      REAL             X( * ), B( * )
      CHARACTER        WHICH*4
*     ..
*     .. Local Scalars ..
      LOGICAL          LSAME
*     ..
*     .. Common Blocks ..
      CHARACTER           CURPFORM*5
      COMMON            / FORMS  / CURPFORM
*
*     .. Executable Statements ..
*
      IF ( LSAME( CURPFORM, 'IDENT' ) ) THEN
         CALL PSOLVENONE( X, B )
      ELSE IF ( LSAME( CURPFORM, 'JACBI' ) ) THEN
         IF( LSAME( WHICH, 'LEFT' ) ) THEN
            CALL PSOLVEJAC( X, B )
         ELSE
            CALL PSOLVENONE( X, B )
         ENDIF
      ELSE
         WRITE(*,*) 'IN PSOLVEQ: UNKNOWN PRECONDITIONER', CURPFORM, 
     $              ' QUITTING' 
         STOP
      ENDIF
*
      RETURN
*
      END
*
      SUBROUTINE QMR( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                MATVECTRANS, PSOLVEQ, PSOLVETRANSQ, INFO )
*
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ
*
*  Purpose
*  =======
*
*  QMR Method solves the linear system Ax = b using the
*  Quasi-Minimal Residual iterative method with preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X      (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*
*  WORK    (workspace) REAL             array, dimension (LDW,11).
*          Workspace for residual, direction vector, etc.
*          Note that W and WTLD, Y and YTLD, and Z and ZTLD share
*          workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  MATVECTRANS  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A'*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A' is the tranpose of a matrix A. Vector x must remain
*          unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVECTRANS( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVEQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. As QMR uses left
*          and right preconditioning and the preconditioners are in
*          common, we must specify in the call which to use. Vector b
*          must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVEQ( X, B, 'LEFT' )
*
*          The preconditioner is passed into the routine in a common block.
*
*  PSOLVETRANSQ  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M'*x = b,
*
*          where x and y are vectors, and M' is the tranpose of a 
*          matrix M. As QMR uses left and right preconditioning and
*          the preconditioners are in common, we must specify in the 
*          call which to use. Vector b must remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVETRANSQ( X, B, 'LEFT' )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO   < BREAKTOL: RHO and RTLD have become
*                                         orthogonal.
*                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
*                                         Convergence has stalled.
*                  -12: GAMMA < BREAKTOL: THETA too large. 
*                                         Convergence has stalled.
*                  -13: DELTA < BREAKTOL: Y and Z have become
*                                         orthogonal.
*                  -14: EPS   < BREAKTOL: Q and PTLD have become
*                                         orthogonal.
*                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  ==============================================================
*
*This variable used to communicate requests between QMR() and 
* QMRREVCOM()
*QMR -> QMRREVCOM: 1 = init, 
*                  2 = use saved state to resume flow.
*QMRREVCOM -> QMR: -1 = done, return to main, 
*                   1 = matvec using SCLR1/2, NDX1/2 
*                   2 = transpose-matvec using SCLR1/2, NDX1/2 
*                   3 = left prec solve using NDX1/2
*                   4 = right prec solve using NDX1/2
*                   5 = left prec transpose-solve using NDX1/2
*                   6 = right transpose-solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*     Arg/Result indices into WORK[].
      INTEGER          NDX1, NDX2
*     Scalars passed from QMRREVCOM to QMR.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         QMRREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == D; 3 == P; 4 == PTLD; 5 == Q; 6 == S; 7 == V;
*     8 == VTLD; 9 == W; 10 == WTLD; 11 == Y; 12 == YTLD; 13 == Z;
*     14 == ZTLD; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL QMRREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM we use the table (BiCGREVCOM -> BiCG)
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec.
            CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call transpose-matvec.
            CALL MATVECTRANS(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 3) THEN
*           call left prec solve.
            CALL PSOLVEQ(WORK(NDX1), WORK(NDX2), 'LEFT')
          ELSEIF (IJOB .eq. 4) THEN
*           call right prec solve.
            CALL PSOLVEQ(WORK(NDX1), WORK(NDX2), 'RIGHT')
          ELSEIF (IJOB .eq. 5) THEN
*           call left prec transpose-solve.
            CALL PSOLVETRANSQ(WORK(NDX1), WORK(NDX2),'LEFT')
          ELSEIF (IJOB .eq. 6) THEN
*           call right prec transpose-solve.
            CALL PSOLVETRANSQ(WORK(NDX1), WORK(NDX2),'RIGHT')
          ELSEIF (IJOB .eq. 7) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
         ELSEIF (IJOB .EQ. 8) THEN
*           do stopping test 2
*           if first time, set INFO so that BNRM2 is computed.
            IF( FTFLG ) INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
            FTFLG = .FALSE.
         ENDIF
*
*        done what revcom asked, set IJOB & go back to it.
         IJOB = 2
         GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
*
      RETURN
*
*     End of QMR
*
      END
*
      SUBROUTINE QMRREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               X( * ), B( * ), WORK( LDW,* )
*     ..
*  Purpose
*  =======
*
*  QMR Method solves the linear system Ax = b using the
*  Quasi-Minimal Residual iterative method with preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X      (input/output) REAL             array, dimension N.
*          On input, the initial guess; on exit, the iterated solution.
*
*
*  WORK    (workspace) REAL             array, dimension (LDW,11).
*          Workspace for residual, direction vector, etc.
*          Note that W and WTLD, Y and YTLD, and Z and ZTLD share
*          workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL            
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO   < BREAKTOL: RHO and RTLD have become
*                                         orthogonal.
*                  -11: BETA  < BREAKTOL: EPS too small in relation to DELTA.
*                                         Convergence has stalled.
*                  -12: GAMMA < BREAKTOL: THETA too large. 
*                                         Convergence has stalled.
*                  -13: DELTA < BREAKTOL: Y and Z have become
*                                         orthogonal.
*                  -14: EPS   < BREAKTOL: Q and PTLD have become
*                                         orthogonal.
*                  -15: XI    < BREAKTOL: Z too small. Convergence has stalled.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*  ==============================================================
*
*     .. Parameters ..
      REAL             ONE, ZERO
      PARAMETER      ( ONE = 1.0E+0 , ZERO = 0.0E+0)
*
*     .. Local Scalars ..
      INTEGER          R, D, P, PTLD, Q, S, V, VTLD, W, WTLD, Y, YTLD,
     $                 Z, ZTLD, MAXIT, NEED1, NEED2
      REAL             TOL, BETA, GAMMA, GAMMA1, DELTA, EPS, ETA, XI, 
     $                 RHO, RHO1, THETA, THETA1, BNRM2,
     $                 RHOTOL, BETATOL, GAMMATOL, DELTATOL, 
     $                 EPSTOL, XITOL, GETBREAK,
     $                 SDOT, SNRM2, C1
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     ..
*     .. External Routines ..
      EXTERNAL         SAXPY, SCOPY, SDOT, SNRM2, SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC        ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
         IF (RLBL .eq. 8) GOTO 8
         IF (RLBL .eq. 9) GOTO 9
         IF (RLBL .eq. 10) GOTO 10
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R     = 1
      D     = 2
      P     = 3
      PTLD  = 4
      Q     = 5
      S     = 6
      V     = 7
      VTLD  = 8
      W     = 9
      WTLD  = 9
      Y     = 10
      YTLD  = 10
      Z     = 11
      ZTLD  = 11
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((D - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((VTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.9 ) THEN
            NEED1 = ((W - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.10 ) THEN
            NEED1 = ((WTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.11 ) THEN
            NEED1 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.12 ) THEN
            NEED1 = ((YTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.13 ) THEN
            NEED1 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.14 ) THEN
            NEED1 = ((ZTLD - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((D - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((PTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((Q - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((VTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.9 ) THEN
            NEED2 = ((W - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.10 ) THEN
            NEED2 = ((WTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.11 ) THEN
            NEED2 = ((Y - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.12 ) THEN
            NEED2 = ((YTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.13 ) THEN
            NEED2 = ((Z - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.14 ) THEN
            NEED2 = ((ZTLD - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set breakdown tolerances.
*
      RHOTOL   = GETBREAK()
      BETATOL  = GETBREAK()
      GAMMATOL = GETBREAK()
      DELTATOL = GETBREAK()
      EPSTOL   = GETBREAK()
      XITOL    = GETBREAK()
*
*     Set initial residual.
*
      CALL SCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ZERO, WORK(1,R) )
*        Note: using D as temp
*********CALL SCOPY( N, X, 1, WORK(1,D), 1 )
         SCLR1 = -ONE
         SCLR2 = ZERO
         NDX1 = ((D - 1) * LDW) + 1
         NDX2 = ((R - 1) * LDW) + 1
         RLBL = 2
         IJOB = 7
         RETURN
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,R), 1 ) .LT. TOL ) GO TO 30
      ENDIF
*
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,VTLD), 1 )
******CALL PSOLVEQ( WORK(1,Y), WORK(1,VTLD), 'LEFT' )
*
         NDX1 = ((Y    - 1) * LDW) + 1
         NDX2 = ((VTLD - 1) * LDW) + 1
         RLBL = 3
         IJOB = 3
         RETURN
*****************
 3       CONTINUE
*****************
*
      RHO   = SNRM2( N, WORK(1,Y), 1 )
*
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,WTLD), 1 )
******CALL PSOLVETRANSQ( WORK(1,Z), WORK(1,WTLD), 'RIGHT' )
*
         NDX1 = ((Z    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 4
         IJOB = 6
         RETURN
*****************
 4       CONTINUE
*****************
*
      XI   = SNRM2( N, WORK(1,Z), 1 )
*
      GAMMA = ONE
      ETA   = -ONE
      THETA = ZERO
*
      ITER = 0
*
   40 CONTINUE
*
*     Perform Preconditioned QMR iteration.
*
         ITER = ITER + 1
*
         IF ( ( ABS( RHO ).LT.RHOTOL ).OR.( ABS( XI ).LT.XITOL ) )
     $      GO TO 25
*
         CALL SCOPY( N, WORK(1,VTLD), 1, WORK(1,V), 1 )
         CALL SSCAL( N, ONE / RHO, WORK(1,V), 1 )
         CALL SSCAL( N, ONE / RHO, WORK(1,Y), 1 )
*
         CALL SCOPY( N, WORK(1,WTLD), 1, WORK(1,W), 1 )
         CALL SSCAL( N, ONE / XI, WORK(1,W), 1 )
         CALL SSCAL( N, ONE / XI, WORK(1,Z), 1 )
*
         DELTA = SDOT( N, WORK(1,Z), 1, WORK(1,Y), 1 )
         IF ( ABS( DELTA ).LT.DELTATOL ) GO TO 25
*
*********CALL PSOLVEQ( WORK(1,YTLD), WORK(1,Y), 'RIGHT' )
*
         NDX1 = ((YTLD - 1) * LDW) + 1
         NDX2 = ((Y    - 1) * LDW) + 1
         RLBL = 5
         IJOB = 4
         RETURN
*****************
 5       CONTINUE
*****************
*
*********CALL PSOLVETRANSQ( WORK(1,ZTLD), WORK(1,Z), 'LEFT' )
*
         NDX1 = ((ZTLD - 1) * LDW) + 1
         NDX2 = ((Z    - 1) * LDW) + 1
         RLBL = 6
         IJOB = 5
         RETURN
*****************
 6       CONTINUE
*****************
*
*
         IF ( ITER.GT.1 ) THEN
            C1 = -( XI * DELTA / EPS )
            CALL SAXPY( N, C1, WORK(1,P), 1, WORK(1,YTLD), 1 )
            CALL SCOPY( N, WORK(1,YTLD), 1, WORK(1,P), 1 )
            CALL SAXPY( N, -( RHO * DELTA / EPS ), WORK(1,Q), 1, 
     $                  WORK(1,ZTLD), 1 )
            CALL SCOPY( N, WORK(1,ZTLD), 1, WORK(1,Q), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,YTLD), 1, WORK(1,P), 1 )
            CALL SCOPY( N, WORK(1,ZTLD), 1, WORK(1,Q), 1 )
         ENDIF
*
*********CALL MATVEC( ONE, WORK(1,P), ZERO, WORK(1,PTLD) )
*
         SCLR1 = ONE
         SCLR2 = ZERO
         NDX1 = ((P    - 1) * LDW) + 1
         NDX2 = ((PTLD - 1) * LDW) + 1
         RLBL = 7
         IJOB = 1
         RETURN
*****************
 7       CONTINUE
*****************
*
*
         EPS = SDOT( N, WORK(1,Q), 1, WORK(1,PTLD), 1 )
         IF ( ABS( EPS ).LT.EPSTOL ) GO TO 25
*
         BETA = EPS / DELTA
         IF ( ABS( BETA ).LT.BETATOL ) GO TO 25
*
         CALL SCOPY( N, WORK(1,PTLD), 1, WORK(1,VTLD), 1 )
         CALL SAXPY( N, -BETA, WORK(1,V), 1, WORK(1,VTLD), 1 )
         CALL PSOLVEQ( WORK(1,Y), WORK(1,VTLD), 'LEFT' )
*
         RHO1 = RHO
         RHO  = SNRM2( N, WORK(1,Y), 1 )
*
         CALL SCOPY( N, WORK(1,W), 1, WORK(1,WTLD), 1 )
*********CALL MATVECTRANS( ONE, WORK(1,Q), -BETA, WORK(1,WTLD) )
*
         SCLR1 = ONE
         SCLR2 = -BETA
         NDX1 = ((Q    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 8
         IJOB = 2
         RETURN
*****************
 8       CONTINUE
*****************
*
*********CALL PSOLVETRANSQ( WORK(1,Z), WORK(1,WTLD), 'RIGHT' )
*
         NDX1 = ((Z    - 1) * LDW) + 1
         NDX2 = ((WTLD - 1) * LDW) + 1
         RLBL = 9
         IJOB = 6
         RETURN
*****************
 9       CONTINUE
*****************
*
*
         XI = SNRM2( N, WORK(1,Z), 1 )
*
         GAMMA1 = GAMMA
         THETA1 = THETA
*
         THETA  = RHO / ( GAMMA1 * ABS( BETA ) )
         GAMMA = ONE / SQRT( ONE + THETA**2 )
         IF ( ABS( GAMMA ).LT.GAMMATOL ) GO TO 25
*
         ETA = -ETA * RHO1 * GAMMA**2 / ( BETA * GAMMA1**2 )
*
         IF ( ITER.GT.1 ) THEN
            CALL SSCAL( N, ( THETA1*GAMMA )**2, WORK(1,D), 1 )
            CALL SAXPY( N, ETA, WORK(1,P), 1, WORK(1,D), 1 )
            CALL SSCAL( N, ( THETA1 * GAMMA )**2, WORK(1,S), 1 )
            CALL SAXPY( N, ETA, WORK(1,PTLD), 1, WORK(1,S), 1 )
         ELSE
            CALL SCOPY( N, WORK(1,P), 1, WORK(1,D), 1 )
            CALL SSCAL( N, ETA, WORK(1,D), 1 )
            CALL SCOPY( N, WORK(1,PTLD), 1, WORK(1,S), 1 )
            CALL SSCAL( N, ETA, WORK(1,S), 1 )
         ENDIF
*
*        Compute current solution vector x.
*
         CALL SAXPY( N, ONE, WORK(1,D), 1, X, 1 )
*
*        Compute residual vector rk, find norm,
*        then check for tolerance.
*
         CALL SAXPY( N, -ONE, WORK(1,S), 1, WORK(1,R), 1 )
*
*********RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
*********IF ( RESID .LE. TOL  ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 10
         IJOB = 8
         RETURN
*
*****************
 10      CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         GO TO 40
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   25 CONTINUE
*
*     Method breakdown.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( BETA ).LT.BETATOL ) THEN
         INFO = -11
      ELSE IF ( ABS( GAMMA ).LT.GAMMATOL ) THEN
         INFO = -12
      ELSE IF ( ABS( DELTA ).LT.DELTATOL ) THEN
         INFO = -13
      ELSE IF ( ABS( EPS ).LT.EPSTOL ) THEN
         INFO = -14
      ELSE IF ( ABS( XI ).LT.XITOL ) THEN
         INFO = -15
      ENDIF
*
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of QMRREVCOM
*
      END
*
      SUBROUTINE SOR( N, B, X, WORK, LDW, ITER, RESID, MATVEC, 
     $                BACKSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( * )
*     .. 
*     .. Function Arguments ..
      EXTERNAL           MATVEC, BACKSOLVE
*
*  Purpose
*  =======
*
*  SOR solves the linear system Ax = b using the Successive 
*  Over-Relaxation iterative method. 
*  The matrix splitting is formed by copying the strict upper 
*  triangular portion of A onto matrix N, stored in WORK. Matrix M 
*  is the lower triangular portion of A.
*  On exit, matrix A and right hand side b are reset to their 
*  original form.
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (input/workspace) REAL             array, dimension (N*(N+3)).
*          The relaxation parameter, OMEGA, should be input in WORK(1).
*          The amount of workspace can be significantly reduced (to 2*N)
*          by customizing the matrix-vector product and backsolve.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL                     
*          On input, the allowable convergence measure for
*          norm( x - x_1 ) / norm( x ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  BACKSOLVE (external subroutine)
*          The user must provide a subroutine to perform the
*          linear system solve
*
*               x := M*x,
*
*          where x is a vector and M is a lower triangular matrix. 
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL BACKSOLVE( N, M, LDM, X )
*
*          The matrix is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: Relaxation parameter OMEGA not in interval (0,2).
*
*     ==========================================================
*
*
*This variable used to communicate requests between SOR() and SORREVCOM()
*SOR -> SORREVCOM: 1 = init, 
*                2 = use saved state to resume flow.
*SORREVCOM -> SOR: -1 = done, return to main, 
*                 1 = matvec using SCLR1/2, NDX1/2 
*                 2 = solve using NDX1/2
      INTEGER          IJOB
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from SORREVCOM to SOR.
      REAL             SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      REAL             TOL, XNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         SORREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ELSE IF ( ( WORK( 1 ).LE.0 ).OR.( WORK( 1 ).GE.2 ) ) THEN
         INFO = -4
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == X1; 2 == TEMP; 3 == MM; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

          CALL SORREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                   NDX1, NDX2, SCLR1, SCLR2, IJOB)

*         On a return from REVCOM() we use the table
*         to decode IJOB.
          IF (IJOB .eq. -1) THEN
*           revcom wants to terminate, so do it.
            GOTO 2
          ELSEIF (IJOB .eq. 1) THEN
*           call matvec with X.
            CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
          ELSEIF (IJOB .eq. 2) THEN
*           call backsolve with X
            CALL BACKSOLVE(N, WORK(NDX1), LDW, X)
         ELSEIF (IJOB .EQ. 3) THEN
*           do stopping test 2
*           XNRM2 shld be recomputed everytime.
            INFO = -1
            CALL STOPTEST2(N, WORK(NDX1), X, XNRM2, RESID, TOL, INFO)
         ENDIF
*
*     done what revcom asked, set IJOB & go back to it.
      IJOB = 2
      GOTO 1
*
*     come here to terminate
 2    CONTINUE
*
*
      RETURN
*
*     End of SOR
*
      END
*
      SUBROUTINE SORREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                     NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO
      REAL               RESID
      INTEGER            NDX1, NDX2
      REAL               SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      REAL               B( * ), X( * ), WORK( LDW,* )
*     .. 
*
*  Purpose
*  =======
*
*  SOR solves the linear system Ax = b using the Successive 
*  Over-Relaxation iterative method. 
*  The matrix splitting is formed by copying the strict upper 
*  triangular portion of A onto matrix N, stored in WORK. Matrix M 
*  is the lower triangular portion of A.
*  On exit, matrix A and right hand side b are reset to their 
*  original form.
*
*  Relative error measured: norm( X - X_1 ) / norm( X ).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) REAL             array, dimension N.
*          On input, the initial guess. This is commonly set to
*          the zero vector.
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (input/workspace) REAL             array, dimension (N*(N+3)).
*          The relaxation parameter, OMEGA, should be input in WORK(1).
*          The amount of workspace can be significantly reduced (to 2*N)
*          by customizing the matrix-vector product and backsolve.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) REAL                     
*          On input, the allowable convergence measure for
*          norm( x - x_1 ) / norm( x ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*            -5: Erroneous NDX1/NDX2 in INIT call.
*            -6: Erroneous RLBL.
*
*          >  0: Convergence to tolerance not achieved. This will be
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -4: Relaxation parameter OMEGA not in interval (0,2).
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) REAL            .
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS:   SAXPY, SCOPY, SNRM2
*     ==========================================================
*
*     .. External Routines ..
      EXTERNAL           SAXPY, SCOPY, SNRM2
*     ..
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER        ( ONE = 1.0E+0 , ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            MM, X1, TEMP, MAXIT, NEED1, NEED2
      REAL               OMEGA, TOL, BNRM2, SNRM2
*     ..
*
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      X1   = 1
      TEMP = 2
      MM   = 3
*
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((TEMP - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((MM - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((X1 - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((TEMP - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((MM - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set relaxation parameter.
*
      OMEGA = WORK( 1,1 )
      IF ( OMEGA .EQ. ZERO ) OMEGA = ONE
*
*     Compute initial residual for ( convergence criteria ).
*
      CALL SCOPY( N, B, 1, WORK(1,X1), 1 )
      IF ( SNRM2( N, X, 1 ) .NE. ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,X1) )
*
         NDX1 = -1
         NDX2 = ((X1   - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = -ONE
         SCLR2 = ONE
         RLBL = 2
         IJOB = 1
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( SNRM2( N, WORK(1,X1), 1 ).LT.TOL ) GO TO 30
      ENDIF
      BNRM2 = SNRM2( N, B, 1 )
      IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
*
*     Matrix A is set to N. WORK(1:N,1:N) is set to MM.
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','SPLIT')
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform SOR iteration
*
      ITER = ITER + 1
*
*        Save the current approximation to X in X1,
*
         CALL SCOPY( N, X, 1, WORK( 1,X1 ), 1 )
*
*        Apply iteration; result is updated approximation vector x
*
         CALL SCOPY( N, B, 1, WORK( 1,TEMP ), 1 )
*********CALL MATVEC( ONE, X, ONE, WORK( 1,TEMP ) )
*
         NDX1 = -1
         NDX2 = ((TEMP - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ONE
         RLBL = 3
         IJOB = 1
         RETURN
*
*****************
 3       CONTINUE
*****************
*
         CALL SCOPY( N, WORK( 1,TEMP ), 1, X, 1 )
*********CALL BACKSOLVE( N, WORK( 1,MM ), LDW, X )
*
         NDX1 = ((MM   - 1) * LDW) + 1
         NDX2 = -1
*        Prepare for return & return
         RLBL = 4
         IJOB = 2
         RETURN
*
*****************
 4       CONTINUE
*****************
*
*        Compute error and check for acceptable convergence.
*
         CALL SAXPY( N, -ONE, X, 1, WORK( 1,X1 ), 1 )
*
*********RESID = SNRM2( N, WORK( 1,X1 ), 1 ) / 
*****$           SNRM2( N, X, 1 )
*********IF ( RESID.LE.TOL ) GO TO 30
*
         NDX1 = NEED1
         NDX2 = NEED2
*        Prepare for resumption & return
         RLBL = 5
         IJOB = 3
         RETURN
*
*****************
 5       CONTINUE
*****************
         IF( INFO.EQ.1 ) GO TO 30
*
         IF ( ITER.EQ.MAXIT ) THEN
            INFO = 1
            GO TO 20
         ENDIF
*
         GO TO 10
*
   20 CONTINUE
*
*     Iteration fails
*     restore B to original form
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','RESTORE')
*
      RLBL = -1
      IJOB = -1
*
      RETURN
*
 30   CONTINUE
*
*     Iteration successful; restore A and B to original form,
*     compute residual norm, and return
*
      CALL MATSPLIT( OMEGA, B, WORK(1,MM), LDW,'SOR','RESTORE')
*
      INFO = 0
      RLBL = -1
      IJOB = -1
*
      RETURN
*
*     End of SORREVCOM
*
      END
*
      SUBROUTINE STOPTEST2( N, R, B, BNRM2, RESID, TOL, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            N, INFO
      REAL               RESID, TOL, BNRM2
*     ..
*     .. Array Arguments ..
      REAL               R( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  Computes the stopping criterion 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  INFO    (output) INTEGER
*          On exit, 1/0 depending on whether stopping criterion
*          was met or not.
*
*  BNRM2   (input/output) REAL            .
*          On first time entry, will be -1.0.
*          On first time exit will contain norm2(B)
*          On all subsequent entry/exit's unchanged.
*
*  RESID   (output) REAL            .
*          On exit, the computed stopping measure.
*
*  TOL     (input) REAL            .
*          On input, the allowable convergence measure.
*
*  R       (input) REAL             array, dimension N.
*          On entry, the residual.
*          Unchanged on exit.
*
*  B       (input) REAL             array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  BLAS CALLS:   SNRM2
*  ============================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. External Routines ..
      REAL             SNRM2
      EXTERNAL         SNRM2
*     ..
*     .. Executable Statements ..
*
      IF( INFO.EQ.-1 ) THEN
         BNRM2 = SNRM2( N, B, 1 )
         IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
      ENDIF
*
      RESID = SNRM2( N, R, 1 ) / BNRM2
*
      INFO = 0
      IF ( RESID.LE.TOL )
     $     INFO = 1
*
      RETURN
*
*     End of STOPTEST2
*
      END

*
*  This file contains routines used by Jacobi, SOR, and Chebyshev:
*
*  Jacobi/SOR:
*
*     MATSPLIT  calls specific matrix splitting routine
*     JACSPLIT
*     SORSPLIT
*     BACKSOLVE
*
*  Chebyshev:
*
*     GETEIG    computes eigenvalue of iteration matrix. 
*
*     =========================================================== 
      SUBROUTINE MATSPLIT( OMEGA, B, WORK, LDW, METHOD, FLAG )
*
      INTEGER            LDW
      REAL               OMEGA, B( * ), WORK( LDW,* )
      CHARACTER          METHOD*1, FLAG*1
      LOGICAL            LSAME
*     ..
*     .. Parameters ..
*
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
*     .. Common Blocks ..
      INTEGER             N, LDA
      REAL                A, M
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM ),
     $                  / MATDIM / N, LDA
*
      IF ( LSAME( METHOD,'JACOBI' ) ) THEN
         CALL JACSPLIT( N, A, LDA, WORK, LDW, FLAG )
      ELSE IF ( LSAME( METHOD, 'SOR' ) ) THEN
         CALL SORSPLIT( OMEGA, N, A, LDA, B, WORK, LDW, FLAG )
      ELSE 
         WRITE(*,*) 'ERROR: UNKNOW METHOD. QUITTING...'
         STOP
      ENDIF
*
      RETURN
*
      END
*
*     ===========================================================
      SUBROUTINE JACSPLIT( N, A, LDA, WORK, LDW, FLAG )
*
      INTEGER            I, J, N, LDA, LDW
      REAL               ZERO, ONE, A( LDA,* ), WORK( LDW,* )
      CHARACTER          FLAG*1
      LOGICAL            LSAME
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*
      IF ( LSAME( FLAG,'SPLIT' ) ) THEN
         DO 20 I = 1, N
            WORK( I,1 ) = ONE / A( I,I )
            A( I,I ) = ZERO
            DO 10 J = 1, N
               A( I,J ) = -A( I,J )
   10       CONTINUE
   20    CONTINUE 
      ELSE IF( LSAME( FLAG,'RECONSTRUCT' ) ) THEN
         DO 40 I = 1, N
            DO 30 J = 1, N
               A( I,J ) = -A( I,J )
   30       CONTINUE
            A( I,I ) = ONE / WORK( I,1 )
   40    CONTINUE
      ELSE
         WRITE(*,*) 'UNKNOWN SPLITTING OPTION. QUITTING...'
         STOP
      ENDIF
*
      RETURN
*
      END
*
*     ===========================================================
      SUBROUTINE SORSPLIT( OMEGA, N, A, LDA, B, WORK, LDW, FLAG )
*  
      INTEGER            I, J, N, LDA, LDW
      REAL               OMEGA, ZERO, ONE,
     $                   B( * ), A( LDA,* ), WORK( LDW,* )
      CHARACTER          FLAG*3
      LOGICAL            LSAME
      PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*
      IF ( LSAME( FLAG,'SPLIT' ) ) THEN
*
*        Set M.
*
         DO 20 I = 1, N
            WORK( I,I ) = A( I,I )
            DO 10 J = 1, I-1
               WORK( I,J ) = OMEGA * A( I,J )
   10       CONTINUE
   20    CONTINUE
*
*        Set NN and B.
*
*        Temporarily store the matrix A in order to reconstruct 
*        the original matrix. Because the lower triangular portion
*        of A must be zeroed, this is the easiest way to deal with it. 
*        This causes the requirement that WORK be N x (2N+3).
*
         DO 40 I = 1, N
            DO 30 J = 1, N
               WORK( I,J+N+3 ) = A( I,J )
   30       CONTINUE
   40    CONTINUE 
*
         DO 60 I = 1, N
            B( I ) = OMEGA * B( I )
            A( I,I ) = ( ONE-OMEGA ) * A( I,I )
            DO 50 J = I+1, N
               A( I,J ) = -OMEGA * A( I,J )
   50       CONTINUE
   60    CONTINUE
*
         DO 80 I = 2, N
            DO 70 J = 1, I-1
               A( I,J ) = ZERO
   70       CONTINUE
   80    CONTINUE
*
      ELSE IF( LSAME( FLAG,'RECONSTRUCT' ) ) THEN
         DO 100 I = 1, N
            B( I ) = B( I ) / OMEGA
            DO 90 J = 1, N
               A( I,J ) = WORK( I,J+N+3 )
   90       CONTINUE
  100    CONTINUE
*
      ELSE
         WRITE(*,*) 'UNKNOWN SPLITTING OPTION. QUITTING...'
         STOP
      ENDIF

      RETURN
*
      END
*
*     ========================================================
      SUBROUTINE BACKSOLVE( N, A, LDA, X )
*
*     .. Argument Declarations ..
      INTEGER             N, LDA
      REAL                X( * ), A( LDA,* )
*
*     Mask to BLAS routine. X overwritten with inv(A)*X
*
      CALL STRSV( 'LOWER','NOTRANS','NONUNIT', N, A, LDA, X, 1 )
*
      RETURN
*
      END
*
*     ===========================================================
      SUBROUTINE GETEIG( WORK, LDW, EIGMAX, EIGMIN )
*
*     .. Argument Declarations ..
      INTEGER            LDW
      REAL               EIGMAX, EIGMIN, WORK( LDW,* )
*
*     This routine using an LAPACK routine for computing all the
*     eigenvalues of the matrix A. This is for testing purposes only,
*     as this is more expensive than a direct dense solver for the
*     linear system. This is for testing purposes only.
*     ..
*     .. Parameters ..
*
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
*     .. Common Blocks ..
      INTEGER             N, LDA
      REAL                A, M
      CHARACTER           CURPFORM*5
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM ),
     $                  / MATDIM / N, LDA
     $                  / FORMS  / CURPFORM
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO
      REAL               ZERO, MATNORM
      PARAMETER        ( ZERO = 0.0E+0 )
      LOGICAL            LSAMEN
*
*     .. Executable Statements ..
*
*     As the matrix A is overwritten in the following routine, we
*     copy it to temporary workspace.
*
      CALL MATCOPY( N, A, LDA, WORK, LDW )
      IF ( LSAMEN( 3, CURPFORM, 'JACBI' ) ) THEN
         DO 10 I = 1, N
            WORK( I,I ) = WORK( I,I ) / M( I )
   10    CONTINUE
      ENDIF
*
*     Call LAPACK eigenvalue routine.
*
      CALL SSYEV('NO_VEC','UPPER', N, WORK, LDW, WORK( 1,N+1 ), 
     $            WORK( 1,N+2 ), (3*N)-1, INFO )
*
      IF ( INFO.NE.0 ) THEN
         WRITE(*,*) 
         WRITE(*,*) 'WARNING: SSYEV COULD NOT COMPUTE ALL EIGENVALUES.'
         WRITE(*,*) 'SETTING EIGMIN/MAX TO DEFAULT VALUES 0 AND |A|'
         WRITE(*,*) 
         EIGMIN = ZERO
         EIGMAX = MATNORM( N, A, LDA )
         RETURN
      ELSE
         EIGMIN = WORK( 1,N+1 )
         EIGMAX = WORK( N,N+1 )
      ENDIF
*
*     Eigenvalues should be positive.
*
      IF ( EIGMIN.LT.ZERO ) THEN
         WRITE(*,*) 'WARNING: MIN EIGENVALUE < 0: SET TO ZERO'
         EIGMIN = ZERO
      ENDIF
      IF ( EIGMAX.LT.EIGMIN ) THEN
         WRITE(*,*) 'WARNING: MAX EIGENVALUE < MIN: SET TO |A|'
         EIGMAX = MATNORM( N, A, LDA )
      ENDIF
*
      RETURN
*
      END
*
*     ======================================================
      SUBROUTINE MATCOPY( N, A, LDA, WORK, LDW )
*
      INTEGER            N, LDA, LDW, I, J
      REAL               A( LDA,* ), WORK( LDW,* )
*
*     This routine copies one matrix to another. Convenient when
*     array dimensions different.
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            WORK( I,J ) = A( I,J )
   10    CONTINUE
   20 CONTINUE
*
      RETURN
*
      END
*     ======================================================
*

*
*  This file contains the single precision fortran BLAS routines needed 
*  by the Templates. Therefore, if you can link to a BLAS library, this
*  file is not needed.
*
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
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
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end
      subroutine scopy(n,sx,incx,sy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
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
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
      real function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
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
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end
*
      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL               TEMP
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMM .
*
      END
      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SGEMV .
*
      END
*     ==============================================================
*
      SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      REAL               ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of SGER  .
*
      END
      real function snrm2 ( n, sx, incx)
      integer i, incx, ix, j, n, next
      real   sx(1),  cutlo, cuthi, hitest, sum, xmax, zero, one
      data   zero, one /0.0e0, 1.0e0/
c
c     euclidean norm of the n-vector stored in sx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c     modified to correct problem with negative increment, 8/21/90.
c     modified to correct failure to update ix, 1/25/92.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 4.441e-16,  1.304e19 /
c
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
      if(incx.lt.0)i = (-n+1)*incx + 1
      ix = 1
c                                                 begin main loop
   20    go to next,(30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 continue
      ix = j
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
c
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j = ix, n
         if(abs(sx(i)) .ge. hitest) go to 100
         sum = sum + sx(i)**2
         i = i + incx
   95 continue
      snrm2 = sqrt( sum )
      go to 300
c
  200 continue
      ix = ix + 1
      i = i + incx
      if( ix .le. n ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
      end
      subroutine srot (n,sx,incx,sy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),stemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
   30 continue
      return
      end
      subroutine sscal(n,sa,sx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      real sa,sx(1)
      integer i,incx,ix,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      do 10 i = 1,n 
        sx(ix) = sa*sx(ix)
        ix = ix + incx 
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end
      subroutine sswap (n,sx,incx,sy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      real sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
      SUBROUTINE SSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SSYMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars,  A is a symmetric matrix and  B and
*  C are  m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  symmetric matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  symmetric  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  symmetric matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  symmetric matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP1, TEMP2
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Set NROWA as the number of rows of A.
*
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( SIDE, 'L' ) )THEN
*
*        Form  C := alpha*A*B + beta*C.
*
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*B*A + beta*C.
*
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( J, K )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
*
      RETURN
*
*     End of SSYMM .
*
      END
*
************************************************************************
*
      SUBROUTINE SSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      REAL               ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYR2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SSYR2 .
*
      END
      SUBROUTINE SSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  SSYR2K  performs one of the symmetric rank 2k operations
*
*     C := alpha*A*B' + alpha*B*A' + beta*C,
*
*  or
*
*     C := alpha*A'*B + alpha*B'*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
*  matrices in the second case.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
*                                        beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
*                                        beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns  of the  matrices  A and B,  and on  entry  with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrices  A and B.  K must be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  k by n  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDB must be at least  max( 1, n ), otherwise  LDB must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL               TEMP1, TEMP2
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYR2K', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  C := alpha*A*B' + alpha*B*A' + C.
*
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
*
*        Form  C := alpha*A'*B + alpha*B'*A + C.
*
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SSYR2K.
*
      END
      SUBROUTINE STRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL               ALPHA
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  STRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of STRMM .
*
      END
*
************************************************************************
*
      SUBROUTINE STRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  STRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of STRMV .
*
      END
      SUBROUTINE STRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  STRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER        ( ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of STRSV .
*
      END
      SUBROUTINE SSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      REAL               ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      REAL               A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     .. Local Scalars ..
      REAL               TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of SSYMV .
*
      END
      integer function isamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      real sx(1),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
      SUBROUTINE XERBLA ( SRNAME, INFO )
*     ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*6        SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the Level 2 BLAS routines.
*
*  It is called by the Level 2 BLAS routines if an input parameter is
*  invalid.
*
*  Installers should consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Parameters
*  ==========
*
*  SRNAME - CHARACTER*6.
*           On entry, SRNAME specifies the name of the routine which
*           called XERBLA.
*
*  INFO   - INTEGER.
*           On entry, INFO specifies the position of the invalid
*           parameter in the parameter-list of the calling routine.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  Written on 20-July-1986.
*
*     .. Executable Statements ..
*
      WRITE (*,99999) SRNAME, INFO
*
      STOP
*
99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,
     $         ' had an illegal value' )
*
*     End of XERBLA.
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
*
*  The routines in this file are from LAPACK. SSYEV is used by the
*  Chebyshev method for computing the extremal eigenvalues of the
*  matrix. All the other routines in this file support this routine.
*  If you can link to the LAPACK library, this file is not needed.
*
      SUBROUTINE SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) REAL array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for SSYTRD returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, J,
     $                   LLWORK, LOPT
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANSY
      EXTERNAL           LSAME, SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SORGTR, SSCAL, SSTEQR, SSTERF, SSYTRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, 3*N-1 ) ) THEN
         INFO = -8
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 3
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL SSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL SSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
      END IF
*
*     Call SSYTRD to reduce symmetric matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL SSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 2*N + WORK( INDWRK )
*
*     For eigenvalues only, call SSTERF.  For eigenvectors, first call
*     SORGTR to generate the orthogonal matrix, then call SSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL SORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL SSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = MAX( 3*N-1, LOPT )
*
      RETURN
*
*     End of SSYEV
*
      END
      SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band symmetric matrix can also be found
*  if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to
*  tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  symmetric matrix.  On entry, Z must contain the
*                  orthogonal matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) REAL array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the orthogonal
*          matrix used in the reduction to tridiagonal form.
*          On exit, if  COMPZ = 'V', Z contains the orthonormal
*          eigenvectors of the original symmetric matrix, and if
*          COMPZ = 'I', Z contains the orthonormal eigenvectors of
*          the symmetric tridiagonal matrix.  If an error exit is
*          made, Z contains the eigenvectors associated with the
*          stored eigenvalues.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) REAL array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is orthogonally similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0, ONE = 1.0, TWO = 2.0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, J, JTOT, K, L, L1, LEND, LENDM1,
     $                   LENDP1, LM1, M, MM, MM1, NM1, NMAXIT
      REAL               B, C, EPS, F, G, P, R, RT1, RT2, S, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLAPY2
      EXTERNAL           LSAME, SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAE2, SLAEV2, SLARTG, SLASR, SLAZRO, SSWAP,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.GT.0 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Determine the unit roundoff for this environment.
*
      EPS = SLAMCH( 'E' )
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL SLAZRO( N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $         GO TO 30
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LEND = M
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         L = LEND
         LEND = L1
      END IF
      L1 = M + 1
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $            GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL SLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL SLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 10
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL SLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 10
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M-1 ) ) ) )
     $            GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL SLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL SLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 10
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL SLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 10
*
      END IF
*
*     Set error -- no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
  140 CONTINUE
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      RETURN
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      DO 180 II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO 170 J = II, N
            IF( D( J ).LT.P ) THEN
               K = J
               P = D( J )
            END IF
  170    CONTINUE
         IF( K.NE.I ) THEN
            D( K ) = D( I )
            D( I ) = P
            IF( ICOMPZ.GT.0 )
     $         CALL SSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
         END IF
  180 CONTINUE
*
      RETURN
*
*     End of SSTEQR
*
      END
      SUBROUTINE SLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  SLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a faster version of the BLAS1 routine SROTG, except for
*  the following differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in SBDSQR when
*        there are zeros on the diagonal).
*
*  Arguments
*  =========
*
*  F       (input) REAL
*          The first component of vector to be rotated.
*
*  G       (input) REAL
*          The second component of vector to be rotated.
*
*  CS      (output) REAL
*          The cosine of the rotation.
*
*  SN      (output) REAL
*          The sine of the rotation.
*
*  R       (output) REAL
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               T, TT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF( ABS( F ).GT.ABS( G ) ) THEN
            T = G / F
            TT = SQRT( ONE+T*T )
            CS = ONE / TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F / G
            TT = SQRT( ONE+T*T )
            SN = ONE / TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
*
*     End of SLARTG
*
      END
      SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), S( * )
*     ..
*
*  Purpose
*  =======
*
*  SLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n real matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  This version vectorises across rows of the array A when SIDE = 'L'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) REAL arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               CTEMP, STEMP, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of SLASR
*
      END
      SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  Purpose
*  =======
*
*  SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) REAL
*          The (1,1) entry of the 2-by-2 matrix.
*
*  B       (input) REAL
*          The (1,2) entry and the conjugate of the (2,1) entry of the
*          2-by-2 matrix.
*
*  C       (input) REAL
*          The (2,2) entry of the 2-by-2 matrix.
*
*  RT1     (output) REAL
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) REAL
*          The eigenvalue of smaller absolute value.
*
*  CS1     (output) REAL
*  SN1     (output) REAL
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      REAL               AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     Compute the eigenvector
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     End of SLAEV2
*
      END
      SUBROUTINE SLAZRO( M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N
      REAL               ALPHA, BETA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SLAZRO initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) REAL
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) REAL
*          The constant to which the diagonal elements are to be set.
*
*  A       (output) REAL array, dimension (LDA,N)
*          On exit, the leading m by n submatrix of A is set such that
*             A(i,j) = ALPHA,  1 <= i <= m, 1 <= j <= n, i <> j
*             A(i,i) = BETA,   1 <= i <= min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ALPHA
   10    CONTINUE
   20 CONTINUE
*
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = BETA
   30 CONTINUE
*
      RETURN
*
*     End of SLAZRO
*
      END
      SUBROUTINE SORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SORGTR generates a real orthogonal matrix Q which is defined as the
*  product of n-1 elementary reflectors of order N, as returned by
*  SSYTRD:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangle of A contains elementary reflectors
*                 from SSYTRD;
*          = 'L': Lower triangle of A contains elementary reflectors
*                 from SSYTRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by SSYTRD.
*          On exit, the N-by-N orthogonal matrix Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  TAU     (input) REAL array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by SSYTRD.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N-1).
*          For optimum performance LWORK >= (N-1)*NB, where NB is
*          the optimal blocksize.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SORGQL, SORGQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to SSYTRD with UPLO = 'U'
*
*        Shift the vectors which define the elementary reflectors one
*        column to the left, and set the last row and column of Q to
*        those of the unit matrix
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL SORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to SSYTRD with UPLO = 'L'.
*
*        Shift the vectors which define the elementary reflectors one
*        column to the right, and set the first row and column of Q to
*        those of the unit matrix
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL SORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      RETURN
*
*     End of SORGTR
*
      END
      SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by SGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by SGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by SGEQRF.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size.
*
      NB = ILAENV( 1, 'SORGQR', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'SORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL SORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL SLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL SORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of SORGQR
*
      END
      SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SORG2R generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by SGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by SGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by SGEQRF.
*
*  WORK    (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARF, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL SSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of SORG2R
*
      END
      SUBROUTINE SORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  SORGQL generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the last N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by SGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by SGEQLF in the last k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by SGEQLF.
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, NB, NBMIN,
     $                   NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORG2L, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORGQL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size.
*
      NB = ILAENV( 1, 'SORGQL', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'SORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL SORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL SLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL SLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL SORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of SORGQL
*
      END
      SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  SLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) REAL array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) REAL array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) REAL array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) REAL array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, STRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of SLARFB
*
      END
      SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      REAL               T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  SLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) REAL array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) REAL array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) REAL array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, STRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL SGEMV( 'Transpose', N-I+1, I-1, -TAU( I ),
     $                        V( I, 1 ), LDV, V( I, I ), 1, ZERO,
     $                        T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  CALL SGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL SGEMV( 'Transpose', N-K+I, K-I, -TAU( I ),
     $                           V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO,
     $                           T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL SGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of SLARFT
*
      END
      SUBROUTINE SORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SORG2L generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the last n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by SGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by SGEQLF in the last k columns of its array
*          argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) REAL array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by SGEQLF.
*
*  WORK    (workspace) REAL array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARF, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL SLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL SSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of SORG2L
*
      END
      SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL               TAU
*     ..
*     .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) REAL array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) REAL
*          The value tau in the representation of H.
*
*  C       (input/output) REAL array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) REAL array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL SGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO,
     $                  WORK, 1 )
*
*           C := C - v * w'
*
            CALL SGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL SGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL SGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of SLARF
*
      END
      SUBROUTINE SSTERF( N, D, E, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      REAL               D( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  SSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) REAL array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) REAL array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm failed to find all of the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0, ONE = 1.0, TWO = 2.0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, JTOT, K, L, L1, LEND, LENDM1, LENDP1,
     $                   LM1, M, MM1, NM1, NMAXIT
      REAL               ALPHA, BB, C, EPS, GAMMA, OLDC, OLDGAM, P, R,
     $                   RT1, RT2, RTE, S, SIGMA, TST
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLAPY2
      EXTERNAL           SLAMCH, SLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAE2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'SSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the unit roundoff for this environment.
*
      EPS = SLAMCH( 'E' )
*
*     Compute the eigenvalues of the tridiagonal matrix.
*
      DO 10 I = 1, N - 1
         E( I ) = E( I )**2
   10 CONTINUE
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   20 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 30 M = L1, NM1
            TST = SQRT( ABS( E( M ) ) )
            IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $         GO TO 40
   30    CONTINUE
      END IF
      M = N
*
   40 CONTINUE
      L = L1
      LEND = M
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         L = LEND
         LEND = L1
      END IF
      L1 = M + 1
*
      IF( LEND.GE.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 60 M = L, LENDM1
               TST = SQRT( ABS( E( M ) ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
*
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        If remaining matrix is 2 by 2, use SLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL SLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 20
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = SLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         MM1 = M - 1
         DO 80 I = MM1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        Eigenvalue found.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 20
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
  100    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 110 M = L, LENDP1, -1
               TST = SQRT( ABS( E( M-1 ) ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M-1 ) ) ) )
     $            GO TO 120
  110       CONTINUE
         END IF
*
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        If remaining matrix is 2 by 2, use SLAE2 to compute its
*        eigenvalues.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL SLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 20
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        Form shift.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = SLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        Inner loop
*
         LM1 = L - 1
         DO 130 I = M, LM1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( LM1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        Eigenvalue found.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 20
*
      END IF
*
*     Set error -- no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
  150 CONTINUE
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      RETURN
*
*     Sort eigenvalues in increasing order.
*
  170 CONTINUE
      DO 190 II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO 180 J = II, N
            IF( D( J ).LT.P ) THEN
               K = J
               P = D( J )
            END IF
  180    CONTINUE
         IF( K.NE.I ) THEN
            D( K ) = D( I )
            D( I ) = P
         END IF
  190 CONTINUE
*
      RETURN
*
*     End of SSTERF
*
      END
      SUBROUTINE SLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B, C, RT1, RT2
*     ..
*
*  Purpose
*  =======
*
*  SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) REAL
*          The (1,1) entry of the 2-by-2 matrix.
*
*  B       (input) REAL
*          The (1,2) and (2,1) entries of the 2-by-2 matrix.
*
*  C       (input) REAL
*          The (2,2) entry of the 2-by-2 matrix.
*
*  RT1     (output) REAL
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) REAL
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
*     ..
*     .. Local Scalars ..
      REAL               AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
*     Compute the eigenvalues
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        Includes case AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        Order of execution important.
*        To get fully accurate smaller eigenvalue,
*        next line needs to be executed in higher precision.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        Includes case RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     End of SLAE2
*
      END
      SUBROUTINE SSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993 
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYTRD reduces a real symmetric matrix A to real symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) REAL array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) REAL array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) REAL array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLATRD, SSYR2K, SSYTD2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size.
*
      NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'SSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = LWORK / LDWORK
               NBMIN = ILAENV( 2, 'SSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL SLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W' - W*V'
*
            CALL SSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL SSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL SLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
*           an update of the form:  A := A - V*W' - W*V'
*
            CALL SSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL SSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of SSYTRD
*
      END
      SUBROUTINE SSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), D( * ), E( * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  SSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
*  form T by an orthogonal similarity transformation: Q' * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the orthogonal
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the orthogonal matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) REAL array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) REAL array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) REAL array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0, ZERO = 0.0, HALF = 1.0 / 2.0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      REAL               ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SLARFG, SSYMV, SSYR2, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            CALL SLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL SSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*SDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL SAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL SSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            CALL SLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL SSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*SDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL SAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL SSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of SSYTD2
*
      END
      SUBROUTINE SLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  SLATRD reduces NB rows and columns of a real symmetric matrix A to
*  symmetric tridiagonal form by an orthogonal similarity
*  transformation Q' * A * Q, and returns the matrices V and W which are
*  needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', SLATRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', SLATRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by SSYTRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements above the diagonal
*            with the array TAU, represent the orthogonal matrix Q as a
*            product of elementary reflectors;
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            with the array TAU, represent the  orthogonal matrix Q as a
*            product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= (1,N).
*
*  E       (output) REAL array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) REAL array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) REAL array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*  and tau in TAU(i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*  and tau in TAU(i).
*
*  The elements of the vectors v together form the n-by-nb matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a symmetric rank-2k update of the form:
*  A := A - V*W' - W*V'.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  a   a   a   v4  v5 )              (  d                  )
*    (      a   a   v4  v5 )              (  1   d              )
*    (          a   1   v5 )              (  v1  1   a          )
*    (              d   1  )              (  v1  v2  a   a      )
*    (                  d  )              (  v1  v2  a   a   a  )
*
*  where d denotes a diagonal element of the reduced matrix, a denotes
*  an element of the original matrix that is unchanged, and vi denotes
*  an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, HALF = 0.5E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      REAL               ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SGEMV, SLARFG, SSCAL, SSYMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               CALL SGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL SGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               CALL SLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL SSYMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL SGEMV( 'Transpose', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL SGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL SGEMV( 'Transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL SGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL SSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*SDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL SAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            CALL SGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL SGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               CALL SLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL SSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL SGEMV( 'Transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL SGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL SGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL SGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL SSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*SDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL SAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of SLATRD
*
      END
      SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               ALPHA, TAU
*     ..
*     .. Array Arguments ..
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  SLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) REAL
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) REAL array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX <> 0.
*
*  TAU     (output) REAL
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2
      EXTERNAL           SLAMCH, SLAPY2, SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = SNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = SLAMCH( 'S' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = SNRM2( N-1, X, INCX )
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of SLARFG
*
      END
      REAL             FUNCTION SLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  SLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by SLAMCH:
*          = 'E' or 'e',   SLAMCH := eps
*          = 'S' or 's ,   SLAMCH := sfmin
*          = 'B' or 'b',   SLAMCH := base
*          = 'P' or 'p',   SLAMCH := eps*base
*          = 'N' or 'n',   SLAMCH := t
*          = 'R' or 'r',   SLAMCH := rnd
*          = 'M' or 'm',   SLAMCH := emin
*          = 'U' or 'u',   SLAMCH := rmin
*          = 'L' or 'l',   SLAMCH := emax
*          = 'O' or 'o',   SLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      SLAMCH = RMACH
      RETURN
*
*     End of SLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  SLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  SLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = SLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of SLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  SLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) REAL
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) REAL
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) REAL
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  SLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine SLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call SLAMC5 to compute EMAX and RMAX.
*
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of SLAMC2
*
      END
*
************************************************************************
*
      REAL             FUNCTION SLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               A, B
*     ..
*
*  Purpose
*  =======
*
*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) REAL
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      SLAMC3 = A + B
*
      RETURN
*
*     End of SLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
*     ..
*
*  Purpose
*  =======
*
*  SLAMC4 is a service routine for SLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) REAL
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of SLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
*     ..
*
*  Purpose
*  =======
*
*  SLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) REAL
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of SLAMC5
*
      END
      REAL             FUNCTION SLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               X, Y
*     ..
*
*  Purpose
*  =======
*
*  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) REAL
*  Y       (input) REAL
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         SLAPY2 = W
      ELSE
         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of SLAPY2
*
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 20, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      REAL             FUNCTION SLANSY( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  SLANSY  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real symmetric matrix A.
*
*  Description
*  ===========
*
*  SLANSY returns the value
*
*     SLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in SLANSY as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is to be referenced.
*          = 'U':  Upper triangular part of A is referenced
*          = 'L':  Lower triangular part of A is referenced
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, SLANSY is
*          set to zero.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The symmetric matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of A contains the upper triangular part
*          of the matrix A, and the strictly lower triangular part of A
*          is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of A contains the lower triangular part of
*          the matrix A, and the strictly upper triangular part of A is
*          not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) REAL array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL               ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is symmetric).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL SLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL SLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL SLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      SLANSY = VALUE
      RETURN
*
*     End of SLANSY
*
      END
      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      REAL               X( * )
*     ..
*
*  Purpose
*  =======
*
*  SLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) REAL
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) REAL
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) REAL
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      REAL               ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of SLASSQ
*
      END

      LOGICAL          FUNCTION LSAMEN( N, CA, CB )
*
*  -- LAPACK auxiliary routine (version 1.0b) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
*     ..
*
*  Purpose
*  =======
*
*  LSAMEN  tests if the first N letters of CA are the same as the
*  first N letters of CB, regardless of case.
*  LSAMEN returns .TRUE. if CA and CB are equivalent except for case
*  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA)
*  or LEN( CB ) is less than N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of characters in CA and CB to be compared.
*
*  CA      (input) CHARACTER*(*)
*  CB      (input) CHARACTER*(*)
*          CA and CB specify two character strings of length at least N.
*          Only the first N characters of each string will be accessed.
*
*     .. Local Scalars ..
      INTEGER            I
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN
*     ..
*     .. Executable Statements ..
*
      LSAMEN = .FALSE.
      IF( LEN( CA ).LT.N .OR. LEN( CB ).LT.N )
     $   GO TO 20
*
*     Do for each character in the two strings.
*
      DO 10 I = 1, N
*
*        Test if the characters are equal using LSAME.
*
         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) )
     $      GO TO 20
*
   10 CONTINUE
      LSAMEN = .TRUE.
*
   20 CONTINUE
      RETURN
*
*     End of LSAMEN
*
      END
*
      REAL FUNCTION SLARAN( ISEED )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  SLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            IPW2
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      SLARAN = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R*
     $         ( REAL( IT4 ) ) ) ) )
      RETURN
*
*     End of SLARAN
*
      END
*
      SUBROUTINE SNSYCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC,
     $                    MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $                    PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL,
     $                    LTEST, NSYRES, NUMTESTS, NUMSUSP, CRITERR )
*
*     .. Scalar Arguments ..
      INTEGER            LDW, LDX, NUMTESTS, NUMSUSP, CRITERR
      REAL               TOL, SCALEDTOL
      LOGICAL            NSYRES
*     ..
*     .. Array Arguments ..
      REAL               X( LDX,* ), B( * ), X0( * ), WORK( * )
      LOGICAL            LTEST( * )
      CHARACTER*5        PFORM( * )
*     ..
*     .. External Functions ..
      EXTERNAL         MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                 PSOLVEQ, PSOLVETRANSQ, BACKSOLVE
*     ..
*     .. Common Blocks ..
*
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      REAL                A, M
      CHARACTER           CURPFORM*5
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
      COMMON            / FORMS  / CURPFORM
*
*  Purpose
*  =======
*
*  Subroutine to test the performance of the nonstationary template
*  kernels on nonsymmetric matrices.
*
*  Generates, solves, and check accuracy of linear systems.
*
*  Algorithms tested:
*
*     4. BiCG
*     5. CGS
*     6. BiCGSTAB
*     7. GMRESm
*     8. QMR
*
*  Various systems are generated. Each method attempts to solve the
*  system to the input TOL in MAXIT iterations. Each method iterates
*  using various preconditioners.
*
*  The result is compared with the normalized scaled residual
*
*     || b - A*x || / ( ||A||||x||*N*TOL ).
*
*  In order to do this, the solution vectors are stored in  matrix
*  X( LDX,* ). Column j contains the solution vector for algorithm j,
*  j as defined above.
*
*  =================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J, MATFORM, NSYS, MAXIT, RESTRT, NPFORMS,
     $                   IPOINTR, IINDX, IVAL, NX, NY, NZ
      REAL               ANORM, MATNORM,
     $                   ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
      LOGICAL            TSTBICG, TSTCGS, TSTBICGS, TSTGMRES, TSTQMR,
     $                   FLAG, LSAMEN
      CHARACTER          AFORM*4, RHS*4, INITIALGUESS*4
*     ..
*     .. Local Arrays ..
      INTEGER            INFO( 9 ), ITER( 9 )
      REAL               RESID( 9 )
*
*     PDE Coefficient functions.
*
      EXTERNAL           ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
*
*     .. Executable Statements ..
*
*
      NPFORMS  = 2
      NUMTESTS = 0
      NUMSUSP  = 0
      CRITERR  = 0
*
      READ(9,*) NSYS
*
*     Check for quick return.
*
      IF ( NSYS.LT.0 ) THEN
         WRITE(*,*) 'ERROR IN NONSYMMETRIC TESTER: NUMBER OF SYSTEMS TO
     $BE GENERATED IS LESS THAN 0'
         RETURN
      ENDIF
*
      FLAG = .FALSE.
      DO 5 I = 4, 8
         IF ( LTEST( I ) ) FLAG = .TRUE.
    5 CONTINUE
      IF ( .NOT.FLAG ) RETURN
*
      TSTBICG   = LTEST( 4 )
      TSTCGS    = LTEST( 5 )
      TSTBICGS  = LTEST( 6 )
      TSTGMRES  = LTEST( 7 )
      TSTQMR    = LTEST( 8 )
*
   10 CONTINUE
*
      DO 60 MATFORM = 1, NSYS
*
         READ(9,*) AFORM, NX, NY, NZ, RHS, INITIALGUESS
*
         N = NX*NY*NZ
         IPOINTR = 1
         IINDX   = IPOINTR + N+1
         IVAL    = IINDX   + 5*N
*
*        The following matrices are generated using a 5- or 7-point
*        stencil using centered differences on a 1d, 2d, or 3d grid,
*        with Dirichlet boundary conditions.
*
*        The last 7 arguments to this routine are the coefficient 
*        functions for the PDE:
*
*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
*           delx ( d u ) + dely (e u) + delz( f u ) + g u
*
         IF ( LSAMEN( 4, AFORM,'PDE1' ) ) THEN
*
*           u_xx + u_yy + au_x + (a_x/2)u for a = 20exp[3.5(x^2 + y^2 )]
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ZEROFUN,
     $                    HENKFUN, ZEROFUN, ZEROFUN, HENKDFUN )
*
*        The following three PDE are from Yang, "PCG-like methods for
*        nonsymmetric linear systems".
*
         ELSE IF ( LSAMEN( 4, AFORM,'PDE2' ) ) THEN
*
*           u_xx + u_yy + u_zz + 1000u_x
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    THOUSFUN, ZEROFUN, ZEROFUN, ZEROFUN )
         ELSE IF ( LSAMEN( 4, AFORM,'PDE3' ) ) THEN
*
*           u_xx + u_yy + u_zz - 10^5x^2(u_x + u_y + u_z )
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    TEN5X2FUN, TEN5X2FUN, TEN5X2FUN, ZEROFUN )
         ELSE IF ( LSAMEN( 4, AFORM,'PDE4' ) ) THEN
*
*           u_xx + u_yy + u_zz + 1000exp(xyz)( u_x + u_y - u_z )
*
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), ONEFUN, ONEFUN, ONEFUN,
     $                    THOUSXFUN, THOUSXFUN, NEGTHOUSXFUN, ZEROFUN )
         ELSE
            WRITE(*,*) AFORM, 'IS AN UNKNOWM MATRIX TYPE'
            GO TO 60
         ENDIF
*
*        Convert to dense form.
*
         CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR), WORK(IINDX),
     $                    N, A, LDA,'ROW', INFO(1) )
         IF ( INFO(1).NE.0 ) THEN
            WRITE(*,81)  N, AFORM, INFO(1)
            WRITE(10,81) N, AFORM, INFO(1)
            GO TO 60
         ENDIF
*
   15    CONTINUE
*
         CALL VECGEN( RHS, N, A, LDA, B, INFO(2) )
         CALL VECGEN( INITIALGUESS, N, A, LDA, X0, INFO(3) )
*
         IF ( INFO(2).NE.0 ) THEN
            WRITE(*,92)  MATFORM
            WRITE(10,92) MATFORM
            GO TO 10
         ELSE IF ( INFO(3).NE.0 ) THEN
            WRITE(*,93)  MATFORM
            WRITE(10,93) MATFORM
            GO TO 10
         ENDIF
*
   20    CONTINUE
*
         MAXIT = 4*N
         ANORM = MATNORM( N, A, LDA )
*
*        Solve system using the various algorithms, using no
*        preconditioning, then diagonal preconditioning.
*
         DO 50 I = 1, NPFORMS
*
*           Initializations.
*
            DO 30 J = 3, 9
               INFO( J ) = 0
               ITER( J ) = MAXIT
               RESID( J ) = TOL
               CALL SCOPY( N, X0, 1, X(1,J), 1 )
   30       CONTINUE
*
            CALL PRECON( N, A, LDA, PFORM( I ), M, INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,94) PFORM(I)
               WRITE(10,94) PFORM(I)
               GO TO 50
            ENDIF
            CURPFORM = PFORM( I )
*
            IF ( TSTBICG )
     $         CALL BICG( N, B, X(1,4), WORK, LDW, ITER(4), RESID(4),
     $                    MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                    INFO(4) )
*
            IF ( TSTCGS )
     $         CALL CGS( N, B, X(1,5), WORK, LDW, ITER(5), RESID(5),
     $                   MATVEC, PSOLVE, INFO(5) )
*
            IF ( TSTBICGS )
     $         CALL BICGSTAB( N, B, X(1,6), WORK, LDW, ITER(6), 
     $                        RESID(6), MATVEC, PSOLVE, INFO(6) )
*
            IF ( TSTGMRES ) THEN
               RESTRT = N / 2
               IF ( RESTRT.EQ.0 ) RESTRT = N
               CALL GMRES( N, B, X(1,7), RESTRT, WORK, LDW,
     $                     WORK((6+RESTRT)*LDW+1), LDW, ITER(7),
     $                     RESID(7), MATVEC, PSOLVE, INFO(7))
            ENDIF
*
            IF ( TSTQMR )
     $         CALL QMR( N, B, X(1,8), WORK, LDW, ITER(8), RESID(8),
     $                   MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ,
     $                   INFO(8) )
*
*           Check for convergence.
*
            DO 40 J = 4, 8
               IF ( INFO( J ).NE.0 ) NUMSUSP = NUMSUSP + 1
   40       CONTINUE
*
*           Write results to file.
*
            CALL RESULT( N, A, LDA, X, LDX, B, WORK(1), 'NSY', PFORM(I),
     $                   ITER, RESID, TOL, INFO, AFORM, ANORM,
     $                   LTEST, SCALEDTOL, NSYRES, CRITERR )
            NUMTESTS = NUMTESTS + 1
*
   50    CONTINUE
*
   60 CONTINUE
*
   81 FORMAT('WARNING: COULD NOT FORM ORDER ', I4, ' ', A4, 
     $       'MATRIX; INFO=',I2)
   92 FORMAT('ERROR: RHS', A4, ' NOT FORMED')
   93 FORMAT('ERROR: INITIAL GUESS', A4, ' NOT FORMED')
   94 FORMAT('ERROR: PRECONDITIONER', A5, ' NOT FORMED')
*
      RETURN
*
*     -- End of SNSYCHK
*
      END
*
      SUBROUTINE SSPDCHK( X, LDX, B, X0, WORK, LDW, PFORM, MATVEC, 
     $                    MATVECTRANS, PSOLVE, PSOLVETRANS, PSOLVEQ,
     $                    PSOLVETRANSQ, BACKSOLVE, TOL, SCALEDTOL,
     $                    LTEST, SPDRES, NUMTESTS, NUMSUSP, CRITERR )
*
*     .. Scalar Arguments ..
      INTEGER            LDW, LDX, NUMTESTS, NUMSUSP, CRITERR
      REAL               TOL, SCALEDTOL
      LOGICAL            SPDRES
*     ..
*     .. Array Arguments ..
      REAL               X( LDX,* ), B( * ), X0( * ), WORK( * )
      LOGICAL            LTEST( * )
      CHARACTER*5        PFORM( * )
*     ..
*     .. External Functions ..
      EXTERNAL           MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                   PSOLVEQ, PSOLVETRANSQ, BACKSOLVE
*     ..
*     .. Common Blocks ..
*     MAXDIM2 = MAXDIM*MAXDIM.
*
      INTEGER             MAXDIM, MAXDIM2
      PARAMETER         ( MAXDIM = 200, MAXDIM2 = 40000 )
*
      INTEGER             N, LDA
      REAL                A, M
      CHARACTER           CURPFORM*5
*
      COMMON            / SYSTEM / A( MAXDIM2 ), M( MAXDIM )
      COMMON            / MATDIM / N, LDA
      COMMON            / FORMS  / CURPFORM
*
*  Purpose
*  =======
*
*  Subroutine to test the performance of the template kernels
*  on symmetric positivie definite matrices.
*
*  Generates, solves, and checks accuracy of linear systems.
*
*  Algorithms tested:
*
*     1. CG
*     2. Cheby
*     3. SOR
*     4. BiCG
*     5. CGS
*     6. BiCGSTAB
*     7. GMRES
*     8. QMR
*     9. Jacobi ( for diagonally dominant Poisson matrices only )
*
*  Various systems are generated. Each method attempts to solve the 
*  system to the input TOL in MAXIT iterations. Each method iterates 
*  using various preconditioners. 
*
*  The result is compared with the normalized scaled residual 
*
*     || b - A*x || / ( ||A||||x||*N*TOL ). 
*
*  In order to do this, the solution vectors are stored in  matrix 
*  X( LDX,* ). Column j contains the solution vector for algorithm j,
*  j as defined above.
*
*  =================================================================
*
*     .. Parameters ..
*
      REAL               ONE
      PARAMETER        ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, K, NX, NY, NZ, MATFORM, NSYS, NPFORMS,
     $                   MAXIT, RESTRT, IVAL, IINDX, IPOINTR
      REAL               ANORM, MATNORM,
     $                   ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
      LOGICAL            TSTCG, TSTCHEBY, TSTSOR, TSTBICG, TSTCGS,
     $                   TSTBICGS, TSTGMRES, TSTQMR, TSTJACOBI,
     $                   FLAG, LSAMEN
      CHARACTER          AFORM*4, RHS*4, INITIALGUESS*4
*     ..
*     .. Local Arrays ..
      INTEGER            INFO( 9 ), ITER( 9 )
      REAL               RESID( 9 )
*     ..
*     .. External Functions ..
*     PDE Coefficient functions.
      EXTERNAL           ONEFUN, ZEROFUN, NEGONEFUN, HENKFUN, HENKDFUN,
     $                   THOUSFUN, THOUSXFUN, TEN5X2FUN, NEGTHOUSXFUN
*     ..
*     .. Executable Statements ..
*
      NPFORMS  = 2
      NUMTESTS = 0
      NUMSUSP  = 0
      CRITERR  = 0
*
      READ(9,*) NSYS
*
*     Check for quick return.
*
      IF ( NSYS.LT.0 ) THEN
         WRITE(*,*) 'ERROR IN NONSYMMETRIC TESTER: NUMBER OF SYSTEMS TO
     $BE GENERATED IS LESS THAN 0'
         RETURN
      ENDIF
*
      FLAG = .FALSE.
      DO 5 I = 1, 9
         IF ( LTEST( I ) ) FLAG = .TRUE.
    5 CONTINUE
      IF ( .NOT.FLAG ) RETURN
*
      TSTCG     = LTEST( 1 )
      TSTCHEBY  = LTEST( 2 )
      TSTSOR    = LTEST( 3 )
      TSTBICG   = LTEST( 4 )
      TSTCGS    = LTEST( 5 )
      TSTBICGS  = LTEST( 6 )
      TSTGMRES  = LTEST( 7 )
      TSTQMR    = LTEST( 8 )
      TSTJACOBI = LTEST( 9 )
*
   10 CONTINUE
*
      DO 60 MATFORM = 1, NSYS
*
         READ(9,*) AFORM, NX, NY, NZ, RHS, INITIALGUESS
*
*        The following two matrices are generated using a 5- or 7-point 
*        stencil using centered differences on a 1d, 2d, or 3d grid, 
*        with Dirichlet boundary conditions.
*
*        The last 7 arguments to this routine are the coefficient 
*        functions for the PDE:
*
*           delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
*           delx ( d u ) + dely (e u) + delz( f u ) + g u
*
         IF ( LSAMEN( 4, AFORM,'F2SH' ) ) THEN
*
*           -u_xx - u_yy
*
            N = NX*NY*NZ
            IPOINTR = 1
            IINDX   = IPOINTR + N+1
            IVAL    = IINDX   + 5*N
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), NEGONEFUN, NEGONEFUN,
     $                    ZEROFUN, ZEROFUN, ZEROFUN, ZEROFUN, ZEROFUN )
            CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR), 
     $                       WORK(IINDX), N, A, LDA,'ROW', INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,81) N, INFO(1)
               WRITE(10,81) N, INFO(1)
               GO TO 60
            ENDIF
*
         ELSE IF ( LSAMEN( 4, AFORM,'F3SH' ) ) THEN
*
*           -u_xx - u_yy - u_zz
*
            N = NX*NY*NZ
            IPOINTR = 1
            IINDX   = IPOINTR + N+1
            IVAL    = IINDX   + 10*N
            CALL GEN57PT( NX, NY, NZ, WORK(IVAL), WORK(IINDX),
     $                    WORK(IPOINTR), NEGONEFUN, NEGONEFUN,
     $                    NEGONEFUN, ZEROFUN, ZEROFUN, ZEROFUN,
     $                    ZEROFUN )
            CALL COMP2DENSE( WORK(IVAL), WORK(IPOINTR),
     $                       WORK(IINDX), N, A, LDA,'ROW', INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,82) N, INFO(1)
               WRITE(10,82) N, INFO(1)
               GO TO 60
            ENDIF
         ELSE IF ( LSAMEN( 4, AFORM,'WATHEN' ) ) THEN
*
*           Form Wathen matrix. 
*           ( Matrix order will be 3*NX*NY + 2*NX + 2*NY + 1 )
*
            K = 0
            CALL WATHEN( NX, NY, K, N, A, LDA, WORK, LDW, INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,83)  3*NX*NY + 2*NX + 2*NY + 1, INFO(1)
               WRITE(10,83) 3*NX*NY + 2*NX + 2*NY + 1, INFO(1)
               GO TO 60
            ENDIF
         ELSE
            WRITE(*,*) AFORM, 'IS AN UNKNOWM MATRIX TYPE'
            GO TO 60
         ENDIF
*
         MAXIT = 4*N
         ANORM = MATNORM( N, A, LDA )
*
*        Form RHS and initial guess vectors.
*
         CALL VECGEN( RHS, N, A, LDA, B, INFO(2) )
         CALL VECGEN( INITIALGUESS, N, A, LDA, X0, INFO(3) )
*
         IF ( INFO(2).NE.0 ) THEN
            WRITE(*,92)  RHS
            WRITE(10,92) RHS
            GO TO 60
         ELSE IF ( INFO(3).NE.0 ) THEN
            WRITE(*,93)  INITIALGUESS
            WRITE(10,93) INITIALGUESS
            GO TO 60
         ENDIF
*
   20    CONTINUE
*
*        Solve system using the various algorithms, using no 
*        preconditioning, then diagonal preconditioning.
*
         DO 50 I = 1, NPFORMS
*
            DO 30 J = 1, 9
               INFO( J ) = 0
               ITER( J )  = MAXIT
               RESID( J ) = TOL
               CALL SCOPY( N, X0, 1, X(1,J), 1 )
   30       CONTINUE
*
            CALL PRECON( N, A, LDA, PFORM( I ), M, INFO(1) )
            IF ( INFO(1).NE.0 ) THEN
               WRITE(*,94) PFORM(I)
               WRITE(10,94) PFORM(I)
               GO TO 50
            ENDIF
            CURPFORM = PFORM( I )
*
            IF ( TSTCG )
     $         CALL CG( N, B, X(1,1), WORK, LDW, ITER(1), RESID(1),
     $                  MATVEC, PSOLVE, INFO(1) )
*
            IF ( TSTCHEBY )
     $         CALL CHEBY( N, B, X(1,2), WORK, LDW, ITER(2), RESID(2),
     $                     MATVEC, PSOLVE, INFO(2) )
*
            IF  ( TSTSOR ) THEN
*
*              Set OMEGA
*
               IF ( MATFORM.EQ.1 ) THEN
*
*                 Guass-Seidel
*
                  WORK(1) = ONE
               ELSE
                  WORK(1) = 1.2E+0
               ENDIF
               CALL SOR( N, B, X(1,3), WORK, LDW, ITER(3),
     $                   RESID(3), MATVEC, BACKSOLVE, INFO(3) )
            ENDIF
*
            IF ( TSTBICG )
     $         CALL BICG( N, B, X(1,4), WORK, LDW, ITER(4), RESID(4),
     $                    MATVEC, MATVECTRANS, PSOLVE, PSOLVETRANS,
     $                    INFO(4) )
*
            IF ( TSTCGS )
     $         CALL CGS( N, B, X(1,5), WORK, LDW, ITER(5), RESID(5),
     $                   MATVEC, PSOLVE, INFO(5) )
*
            IF ( TSTBICGS )
     $         CALL BICGSTAB( N, B, X(1,6), WORK, LDW, ITER(6),
     $                        RESID(6), MATVEC, PSOLVE, INFO(6) )
*
            IF ( TSTGMRES ) THEN
*
*              For the symmetric case, restarts = N.
*
               RESTRT = N
               CALL GMRES( N, B, X(1,7), RESTRT, WORK, LDW,
     $                     WORK((6+RESTRT)*LDW+1), LDW, ITER(7),
     $                     RESID(7), MATVEC, PSOLVE, INFO(7) )
            ENDIF
*
            IF ( TSTQMR )
     $          CALL QMR( N, B, X(1,8), WORK, LDW, ITER(8), RESID(8),
     $                    MATVEC, MATVECTRANS, PSOLVEQ, PSOLVETRANSQ, 
     $                    INFO(8))
*
            IF ( TSTJACOBI ) THEN
               IF ( .NOT.LSAMEN( 3, AFORM,'WATH' ).AND.( I.EQ.1 ) ) THEN
*
*                 Since preconditioning does not apply to Jacobi, it is
*                 only called once per test matrix. (Wathen matrix is 
*                 not diagonally dominant.)
*
                  CALL JACOBI( N, B, X(1,9), WORK, LDW, ITER(9), 
     $                         RESID(9), MATVEC, INFO(9) )
                  IF ( INFO( 9 ).NE.0 ) NUMSUSP = NUMSUSP + 1
               ELSE
*
*                 Flag not to check accuracy.
*
                  INFO( 9 ) = 100
               ENDIF
            ENDIF
*
*           Check for convergence.
*
            DO 40 J = 1, 8
               IF ( INFO( J ).NE.0 ) NUMSUSP = NUMSUSP + 1 
   40       CONTINUE
*
*           Write results to file.
*
            CALL RESULT( N, A, LDA, X, LDX, B, WORK(1), 'SPD', PFORM(I),
     $                   ITER, RESID, TOL, INFO, AFORM, ANORM,
     $                   LTEST, SCALEDTOL, SPDRES, CRITERR )
            NUMTESTS = NUMTESTS + 1
*
   50    CONTINUE
*
   60 CONTINUE
*
   81 FORMAT('WARNING: COULD NOT FORM 2-D POISSON (N=', I4,'),INFO=',I2)
   82 FORMAT('WARNING: COULD NOT FORM 3-D POISSON (N=', I4,'),INFO=',I2)
   83 FORMAT('WARNING: COULD NOT FORM ORDER', I4,' WATHEN MATRIX')
*
   92 FORMAT('ERROR: RHS FOR TEST MATRIX', A4, ' NOT FORMED')
   93 FORMAT('ERROR: INITIAL GUESS FOR TEST MATRIX', A4, ' NOT FORMED')
   94 FORMAT('ERROR: PRECONDITIONER', A5, ' NOT FORMED')
*
  200 CONTINUE
*
      RETURN
*
*     -- End of SSPDCHK
*
      END
*  This file contains various routines needed by the tester for generating
*  the test matrices. They are all copied from the LAPACK test suite,
*  so if you can link to this tester, this file is not needed.
*
      REAL FUNCTION SLARAN2( ISEED )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  SLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            IPW2
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      SLARAN = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R*
     $         ( REAL( IT4 ) ) ) ) )
      RETURN
*
*     End of SLARAN
*
      END
