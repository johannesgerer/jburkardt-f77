      SUBROUTINE CHRCNT(XMESS,NMESS)
C***BEGIN PROLOGUE  CHRCNT
C***REFER TO  SQED
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CHRCNT
C     COUNT THE NUMBER OF NON-BLANK CHARACTERS IN XMESS.
C     RETURN THE RESULT IN NMESS.
      CHARACTER * (*)  XMESS
      INTEGER NMESS
      DO 10 I = LEN(XMESS),1,-1
         IF (XMESS(I:I).NE.' ') THEN
             NMESS = I
             RETURN

         END IF

   10 CONTINUE
      NMESS = 0
      RETURN
      END
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    D1MACH ( 3) = B**(-T), the smallest relative spacing.
c    D1MACH ( 4) = B**(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer i

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      else if ( i == 1 ) then
        d1mach = 4.450147717014403D-308
      else if ( i == 2 ) then
        d1mach = 8.988465674311579D+307
      else if ( i == 3 ) then
        d1mach = 1.110223024625157D-016
      else if ( i == 4 ) then
        d1mach = 2.220446049250313D-016
      else if ( i == 5 ) then
        d1mach = 0.301029995663981D+000
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      end if

      return
      end
      function dasum ( n, dx, incx )

c*********************************************************************72
c
cc DASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c  Modified:
c
c    08 July 2007
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
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of X.
c    INCX must not be negative.
c
c    Output, double precision DASUM, the sum of the absolute values of X.
c
      implicit none

      double precision dasum
      double precision dtemp
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
      end do
      dasum = dtemp
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dtemp + dabs(dx(i))
      end do
      if( n .lt. 6 ) go to 60

   40 continue

      do i = m + 1, n, 6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     &  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
      end do

   60 dasum = dtemp
      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    08 July 2007
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
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

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
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dy(i) = dy(i) + da*dx(i)
      end do

      if( n .lt. 4 ) return

   40 continue

      do i = m + 1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
      end do

      return
      end
      SUBROUTINE DBOCLS(W,MDW,MCON,MROWS,NCOLS,BL,BU,IND,IOPT,X,RNORMC,
     *   RNORM,MODE,RW,IW)
C***BEGIN PROLOGUE  DBOCLS
C***DATE WRITTEN   821220   (YYMMDD)
C***REVISION DATE  880722   (YYMMDD)
C***CATEGORY NO.  K1A2A,G2E,G2H1,G2H2
C***KEYWORDS  BOUNDS,CONSTRAINTS,INEQUALITY,LEAST SQUARES,LINEAR
C***AUTHOR  HANSON, R. J., SNLA
C***PURPOSE  Solve the bounded and constrained least squares
C            problem consisting of solving the equation
C                      E*X = F  (in the least squares sense)
C             subject to the linear constraints
C                            C*X = Y.
C***DESCRIPTION
C
C     This subprogram solves the bounded and constrained least squares
C     problem. The problem statement is:
C
C     Solve E*X = F (least squares sense), subject to constraints
C     C*X=Y.
C
C     In this formulation both X and Y are unknowns, and both may
C     have bounds on any of their components.  This formulation
C     of the problem allows the user to have equality and inequality
C     constraints as well as simple bounds on the solution components.
C
C     This constrained linear least squares subprogram solves E*X=F
C     subject to C*X=Y, where E is MROWS by NCOLS, C is MCON by NCOLS.
C
C      The user must have dimension statements of the form
C
C      DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON), BU(NCOLS+MCON),
C     * X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
C       INTEGER IND(NCOLS+MCON), IOPT(17+NI), IW(2*(NCOLS+MCON))
C
C     (here NX=number of extra locations required for the options; NX=0
C     if no options are in use. Also NI=number of extra locations
C     for options 1-9.
C
C    INPUT
C    -----
C
C    -------------------------
C    W(MDW,*),MCON,MROWS,NCOLS
C    -------------------------
C     The array W contains the (possibly null) matrix [C:*] followed by
C     [E:F].  This must be placed in W as follows:
C          [C  :  *]
C     W  = [       ]
C          [E  :  F]
C     The (*) after C indicates that this data can be undefined. The
C     matrix [E:F] has MROWS rows and NCOLS+1 columns. The matrix C is
C     placed in the first MCON rows of W(*,*) while [E:F]
C     follows in rows MCON+1 through MCON+MROWS of W(*,*). The vector F
C     is placed in rows MCON+1 through MCON+MROWS, column NCOLS+1. The
C     values of MDW and NCOLS must be positive; the value of MCON must
C     be nonnegative. An exception to this occurs when using option 1
C     for accumulation of blocks of equations. In that case MROWS is an
C     OUTPUT variable only, and the matrix data for [E:F] is placed in
C     W(*,*), one block of rows at a time. See IOPT(*) contents, option
C     number 1, for further details. The row dimension, MDW, of the
C     array W(*,*) must satisfy the inequality:
C
C     If using option 1,
C                     MDW .ge. MCON + max(max. number of
C                     rows accumulated, NCOLS)
C
C     If using option 8, MDW .ge. MCON + MROWS.
C     Else, MDW .ge. MCON + max(MROWS, NCOLS).
C
C     Other values are errors, but this is checked only when using
C     option=2.  The value of MROWS is an output parameter when
C     using option number 1 for accumulating large blocks of least
C     squares equations before solving the problem.
C     See IOPT(*) contents for details about option 1.
C
C    ------------------
C    BL(*),BU(*),IND(*)
C    ------------------
C     These arrays contain the information about the bounds that the
C     solution values are to satisfy. The value of IND(J) tells the
C     type of bound and BL(J) and BU(J) give the explicit values for
C     the respective upper and lower bounds on the unknowns X and Y.
C     The first NVARS entries of IND(*), BL(*) and BU(*) specify
C     bounds on X; the next MCON entries specify bounds on Y.
C
C    1.    For IND(J)=1, require X(J) .ge. BL(J);
C          IF J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J).
C          (the value of BU(J) is not used.)
C    2.    For IND(J)=2, require X(J) .le. BU(J);
C          IF J.gt.NCOLS,        Y(J-NCOLS) .le. BU(J).
C          (the value of BL(J) is not used.)
C    3.    For IND(J)=3, require X(J) .ge. BL(J) and
C                                X(J) .le. BU(J);
C          IF J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J) and
C                                Y(J-NCOLS) .le. BU(J).
C          (to impose equality constraints have BL(J)=BU(J)=
C          constraining value.)
C    4.    For IND(J)=4, no bounds on X(J) or Y(J-NCOLS) are required.
C          (the values of BL(J) and BU(J) are not used.)
C
C     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
C     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J)
C     is  an  error.   The values BL(J), BU(J), J .gt. NCOLS, will be
C     changed.  Significant changes mean that the constraints are
C     infeasible.  (Users must make this decision themselves.)
C     The new values for BL(J), BU(J), J .gt. NCOLS, define a
C     region such that the perturbed problem is feasible.  If users
C     know that their problem is feasible, this step can be skipped
C     by using option number 8 described below.
C
C    -------
C    IOPT(*)
C    -------
C     This is the array where the user can specify nonstandard options
C     for DBOCLS( ). Most of the time this feature can be ignored by
C     setting the input value IOPT(1)=99. Occasionally users may have
C     needs that require use of the following subprogram options. For
C     details about how to use the options see below: IOPT(*) CONTENTS.
C
C     Option Number   Brief Statement of Purpose
C     ------ ------   ----- --------- -- -------
C           1         Return to user for accumulation of blocks
C                     of least squares equations.  The values
C                     of IOPT(*) are changed with this option.
C                     The changes are updates to pointers for
C                     placing the rows of equations into position
C                     for processing.
C           2         Check lengths of all arrays used in the
C                     subprogram.
C           3         Column scaling of the data matrix, [C].
C                                                        [E]
C           4         User provides column scaling for matrix [C].
C                                                             [E]
C           5         Provide option array to the low-level
C                     subprogram DBOLS( ).
C                     {Provide option array to the low-level
C                     subprogram DBOLSM( ) by imbedding an
C                     option array within the option array to
C                     DBOLS(). Option 6 is now disabled.}
C           7         Move the IOPT(*) processing pointer.
C           8         Do not preprocess the constraints to
C                     resolve infeasibilities.
C           9         Do not pretriangularize the least squares matrix.
C          99         No more options to change.
C
C    ----
C    X(*)
C    ----
C     This array is used to pass data associated with options 4,5 and
C     6. Ignore this parameter (on input) if no options are used.
C     Otherwise see below: IOPT(*) CONTENTS.
C
C
C    OUTPUT
C    ------
C
C    -----------------
C    X(*),RNORMC,RNORM
C    -----------------
C     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
C     the constrained least squares problem. The value RNORMC is the
C     minimum residual vector length for the constraints C*X - Y = 0.
C     The value RNORM is the minimum residual vector length for the
C     least squares equations. Normally RNORMC=0, but in the case of
C     inconsistent constraints this value will be nonzero.
C     The values of X are returned in the first NVARS entries of X(*).
C     The values of Y are returned in the last MCON entries of X(*).
C
C    ----
C    MODE
C    ----
C     The sign of MODE determines whether the subprogram has completed
C     normally, or encountered an error condition or abnormal status. A
C     value of MODE .ge. 0 signifies that the subprogram has completed
C     normally. The value of mode (.ge. 0) is the number of variables
C     in an active status: not at a bound nor at the value zero, for
C     the case of free variables. A negative value of MODE will be one
C     of the cases (-57)-(-41), (-37)-(-22), (-19)-(-2). Values .lt. -1
C     correspond to an abnormal completion of the subprogram. These
C     error messages are in groups for the subprograms DBOCLS(),
C     DBOLSM(), and DBOLS().  An approximate solution will be returned
C     to the user only when max. iterations is reached, MODE=-22.
C
C    -----------
C    RW(*),IW(*)
C    -----------
C     These are working arrays.  (normally the user can ignore the
C     contents of these arrays.)
C
C    IOPT(*) CONTENTS
C    ------- --------
C     The option array allows a user to modify some internal variables
C     in the subprogram without recompiling the source code. A central
C     goal of the initial software design was to do a good job for most
C     people. Thus the use of options will be restricted to a select
C     group of users. The processing of the option array proceeds as
C     follows: a pointer, here called LP, is initially set to the value
C     1. At the pointer position the option number is extracted and
C     used for locating other information that allows for options to be
C     changed. The portion of the array IOPT(*) that is used for each
C     option is fixed; the user and the subprogram both know how many
C     locations are needed for each option. The value of LP is updated
C     for each option based on the amount of storage in IOPT(*) that is
C     required. A great deal of error checking is done by the
C     subprogram on the contents of the option array. Nevertheless it
C     is still possible to give the subprogram optional input that is
C     meaningless. For example option 4 uses the locations
C     X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing scaling data.
C     The user must manage the allocation of these locations.
C
C   1
C   -
C     This option allows the user to solve problems with a large number
C     of rows compared to the number of variables. The idea is that the
C     subprogram returns to the user (perhaps many times) and receives
C     new least squares equations from the calling program unit.
C     Eventually the user signals "that's all" and a solution is then
C     computed. The value of MROWS is an output variable when this
C     option is used. Its value is always in the range 0 .le. MROWS
C     .le. NCOLS+1. It is the number of rows after the
C     triangularization of the entire set of equations. If LP is the
C     processing pointer for IOPT(*), the usage for the sequential
C     processing of blocks of equations is
C
C
C        IOPT(LP)=1
C         Move block of equations to W(*,*) starting at
C         the first row of W(*,*).
C        IOPT(LP+3)=# of rows in the block; user defined
C
C     The user now calls DBOCLS( ) in a loop. The value of IOPT(LP+1)
C     directs the user's action. The value of IOPT(LP+2) points to
C     where the subsequent rows are to be placed in W(*,*). Both of
C     these values are first defined in the subprogram. The user
C     changes the value of IOPT(LP+1) (to 2) as a signal that all of
C     the rows have been processed.
C
C
C      .<LOOP
C      . CALL DBOCLS( )
C      . IF(IOPT(LP+1) .EQ. 1) THEN
C      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
C      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
C      .    W(*,*) STARTING AT ROW MCON + IOPT(LP+2).
C      .
C      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
C      .       IOPT(LP+1)=2
C      .<------CYCLE LOOP
C      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN
C      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0
C      . ELSE
C      . ERROR CONDITION; SHOULD NOT HAPPEN.
C      .<END LOOP
C
C     Use of this option adds 4 to the required length of IOPT(*).
C
C   2
C   -
C     This option is useful for checking the lengths of all arrays used
C     by DBOCLS( ) against their actual requirements for this problem.
C     The idea is simple: the user's program unit passes the declared
C     dimension information of the arrays. These values are compared
C     against the problem-dependent needs within the subprogram. If any
C     of the dimensions are too small an error message is printed and a
C     negative value of MODE is returned, -41 to -47. The printed error
C     message tells how long the dimension should be. If LP is the
C     processing pointer for IOPT(*),
C
C        IOPT(LP)=2
C        IOPT(LP+1)=Row dimension of W(*,*)
C        IOPT(LP+2)=Col. dimension of W(*,*)
C        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
C        IOPT(LP+4)=Dimension of X(*)
C        IOPT(LP+5)=Dimension of RW(*)
C        IOPT(LP+6)=Dimension of IW(*)
C        IOPT(LP+7)=Dimension of IOPT(*)
C         .
C        CALL DBOCLS( )
C
C     Use of this option adds 8 to the required length of IOPT(*).
C
C   3
C   -
C     This option can change the type of scaling for the data matrix.
C     Nominally each nonzero column of the matrix is scaled so that the
C     magnitude of its largest entry is equal to the value ONE. If LP
C     is the processing pointer for IOPT(*),
C
C        IOPT(LP)=3
C        IOPT(LP+1)=1,2 or 3
C            1= Nominal scaling as noted;
C            2= Each nonzero column scaled to have length ONE;
C            3= Identity scaling; scaling effectively suppressed.
C         .
C        CALL DBOCLS( )
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   4
C   -
C     This options allows the user to provide arbitrary (positive)
C     column scaling for the matrix. If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=4
C        IOPT(LP+1)=IOFF
C        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
C        = Positive scale factors for cols. of E.
C         .
C        CALL DBOCLS( )
C
C     Use of this option adds 2 to the required length of IOPT(*)
C     and NCOLS to the required length of X(*).
C
C   5
C   -
C     This option allows the user to provide an option array to the
C     low-level subprogram DBOLS( ). If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=5
C        IOPT(LP+1)= Position in IOPT(*) where option array
C                    data for DBOLS( ) begins.
C         .
C        CALL DBOCLS( )
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   6
C   -
C     This option is no longer operative.  To pass an option array
C     to the low-level subprogram DBOLSM( ), imbed it within an option
C     array passed to DBOLS() using option 5.
C
C   7
C   -
C     Move the processing pointer (either forward or backward) to the
C     location IOPT(LP+1). The processing pointer moves to locations
C     LP+2 if option number 7 is used with the value -7.  For
C     example to skip over locations 3,...,NCOLS+2,
C
C       IOPT(1)=7
C       IOPT(2)=NCOLS+3
C       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
C       IOPT(NCOLS+3)=99
C       CALL DBOCLS( )
C
C     CAUTION: Misuse of this option can yield some very hard-to-find
C     bugs. Use it with care. It is intended to be used for passing
C     option arrays to other subprograms.
C
C   8
C   -
C     This option allows the user to suppress the algorithmic feature
C     of DBOCLS( ) that processes the constraint equations C*X = Y and
C     resolves infeasibilities. The steps normally done are to solve
C     C*X - Y = 0 in a least squares sense using the stated bounds on
C     both X and Y. Then the "reachable" vector Y = C*X is computed
C     using the solution X obtained. Finally the stated bounds for Y are
C     enlarged to include C*X. To suppress the feature:
C
C
C       IOPT(LP)=8
C         .
C       CALL DBOCLS( )
C
C     Use of this option adds 1 to the required length of IOPT(*).
C
C   9
C   -
C     This option allows the user to suppress the pretriangularizing
C     step of the least squares matrix that is done within DBOCLS( ).
C     This is primarily a means of enhancing the subprogram efficiency
C     and has little effect on accuracy. To suppress the step, set:
C
C       IOPT(LP)=9
C         .
C       CALL DBOCLS( )
C
C     Use of this option adds 1 to the required length of IOPT(*).
C
C   99
C   --
C     There are no more options to change.
C
C     Only option numbers -99, -9,-8,...,-1, 1,2,...,9, and 99 are
C     permitted. Other values are errors. Options -99,-1,...,-9 mean
C     that the respective options 99,1,...,9 are left at their default
C     values. An example is the option to suppress the preprocessing of
C     contraints:
C
C       IOPT(1)=-8 Option is recognized but not changed
C       IOPT(2)=99
C       CALL DBOCLS( )
C
C    Error Messages for DBOCLS()
C    ----- -------- --- --------
C
C WARNING in...
C DBOCLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE NUMBER
C OF EFFECTIVE ROWS=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        41
C
C WARNING IN...
C DBOCLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+
C MCON+1=(I2).
C           IN ABOVE MESSAGE, I1=         2
C           IN ABOVE MESSAGE, I2=         3
C ERROR NUMBER =        42
C
C WARNING IN...
C DBOCLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1)
C MUST BE .GE. NCOLS+MCON=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        43
C
C WARNING IN...
C DBOCLS(). THE DIMENSION OF X()=(I1) MUST BE
C .GE. THE REQD.LENGTH=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        44
C
C WARNING IN...
C DBOCLS(). THE .
C DBOCLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS+2*MCON=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         4
C ERROR NUMBER =        46
C
C WARNING IN...
C DBOCLS(). THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD.
C LEN.=(I2).
C           IN ABOVE MESSAGE, I1=        16
C           IN ABOVE MESSAGE, I2=        18
C ERROR NUMBER =        47
C
C WARNING IN...
C DBOCLS(). ISCALE OPTION=(I1) MUST BE 1-3.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =        48
C
C WARNING IN...
C DBOCLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED COLUMN SCALING
C MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =        49
C
C WARNING IN...
C DBOCLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE.
C  COMPONENT (I1) NOW = (R1).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, R1=    0.
C ERROR NUMBER =        50
C
C WARNING IN...
C DBOCLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED.
C           IN ABOVE MESSAGE, I1=      1001
C ERROR NUMBER =        51
C
C WARNING IN...
C DBOCLS(). NO. OF ROWS=(I1) MUST BE .GE. 0 .AND. .LE. MDW-MCON=(I2).
C           IN ABOVE MESSAGE, I1=         2
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        52
C
C WARNING IN...
C DBOCLS(). MDW=(I1) MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =        53
C
C WARNING IN...
C DBOCLS(). MCON=(I1) MUST BE NONNEGATIVE.
C           IN ABOVE MESSAGE, I1=        -1
C ERROR NUMBER =        54
C
C WARNING IN...
C DBOCLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =        55
C
C WARNING IN...
C DBOCLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         0
C ERROR NUMBER =        56
C
C WARNING IN...
C DBOCLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, R1=     .1000000000E+01
C           IN ABOVE MESSAGE, R2=    0.
C ERROR NUMBER =        57
C           LINEAR CONSTRAINTS, SNLA REPT. SAND82-1517, AUG., (1982).
C***REFERENCES  HANSON, R. J. LINEAR LEAST SQUARES WITH BOUNDS AND
C                 LINEAR CONSTRAINTS, SIAM J. SCI. STAT. COMPUT., VOL. 7,
C                 NO. 3, JULY, 1986.
C***ROUTINES CALLED  D1MACH,DASUM,DBOLS,DCOPY,DDOT,DNRM2,DSCAL,XERRWV
C***COMMON BLOCKS    DBOCOM
C***END PROLOGUE  DBOCLS
C     REVISED 880722-1100
C     REVISED YYMMDD-HHMM
C
C    PURPOSE
C    -------
C     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE LEAST SQUARES
C     PROBLEM CONSISTING OF LINEAR CONSTRAINTS
C
C              C*X = Y
C
C     AND LEAST SQUARES EQUATIONS
C
C              E*X = F
C
C     IN THIS FORMULATION THE VECTORS X AND Y ARE BOTH UNKNOWNS.
C     FURTHER, X AND Y MAY BOTH HAVE USER-SPECIFIED BOUNDS ON EACH
C     COMPONENT.  THE USER MUST HAVE DIMENSION STATEMENTS OF THE
C     FORM
C
C     DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON),BU(NCOLS+MCON),
C               X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON)
C
C     INTEGER IND(NCOLS+MCON), IOPT(16+NI), IW(2*(NCOLS+MCON))
C
C     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
C     EDITING AT THE CARD 'C++'.
C     CHANGE THIS SUBPROGRAM TO DBOCLS AND THE STRINGS
C     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/, /SRELPR/ TO /DRELPR/,
C     /R1MACH/ TO /D1MACH/, /E0/ TO /D0/, /SCOPY/ TO /DCOPY/,
C     /SSCAL/ TO /DSCAL/, /SASUM/ TO /DASUM/, /DBOLS/ TO /DBOLS/,
C     /REAL            / TO /DOUBLE PRECISION/.
C ++
      DOUBLE PRECISION W(MDW,*),BL(*),BU(*),X(*),RW(*)
      DOUBLE PRECISION ANORM, CNORM, ONE, RNORM, RNORMC, DRELPR
      DOUBLE PRECISION T, T1, T2, DDOT, DNRM2, WT, ZERO
      DOUBLE PRECISION DASUM, D1MACH
C     THIS VARIABLE REMAINS TYPED REAL.
      REAL RDUM
      INTEGER IND(*),IOPT(*),IW(*),JOPT(05)
      LOGICAL CHECKL,FILTER,ACCUM,PRETRI
      COMMON /DBOCOM/ IDOPE(5)
      SAVE IGO,ACCUM,CHECKL
      DATA IGO/0/
C***FIRST EXECUTABLE STATEMENT  DBOCLS
      NERR = 0
      MODE = 0
      LEVEL = 1
      IF (IGO.EQ.0) THEN
C     DO(CHECK VALIDITY OF INPUT DATA)
C     PROCEDURE(CHECK VALIDITY OF INPUT DATA)
C
C     SEE THAT MDW IS .GT.0. GROSS CHECK ONLY.
          IF (MDW.LE.0) THEN
              NERR = 53
              NCHAR = 36
              CALL XERRWV('DBOCLS(). MDW=(I1) MUST BE POSITIVE.',NCHAR,
     *                    NERR,LEVEL,1,MDW,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 260
          END IF
C
C     SEE THAT NUMBER OF CONSTRAINTS IS NONNEGATIVE.
          IF (MCON.LT.0) THEN
              NERR = 54
              NCHAR = 40
              CALL XERRWV('DBOCLS(). MCON=(I1) MUST BE NONNEGATIVE.',
     *                    NCHAR,NERR,LEVEL,1,MCON,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 260
          END IF
C
C     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
          IF (NCOLS.LE.0) THEN
              NERR = 55
              NCHAR = 59
              CALL XERRWV(
     *     'DBOCLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.'
     *                    ,NCHAR,NERR,LEVEL,1,NCOLS,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 260
          END IF
C
C     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
          DO 10 J = 1,NCOLS + MCON
              IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
                  NERR = 56
                  NCHAR = 46
                  CALL XERRWV(
     *                  'DBOCLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.'
     *                        ,NCHAR,NERR,LEVEL,2,J,IND(J),0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
   10     CONTINUE
C
C     SEE THAT BOUNDS ARE CONSISTENT.
          DO 20 J = 1,NCOLS + MCON
              IF (IND(J).EQ.3) THEN
                  IF (BL(J).GT.BU(J)) THEN
                      NERR = 57
                      NCHAR = 58
                      CALL XERRWV(
     *      'DBOCLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2).'
     *                            ,NCHAR,NERR,LEVEL,1,J,IDUM,2,BL(J),
     *                            BU(J))
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 260
                  END IF
              END IF
   20     CONTINUE
C     END PROCEDURE
C     DO(PROCESS OPTION ARRAY)
C     PROCEDURE(PROCESS OPTION ARRAY)
          ZERO = 0.D0
          ONE = 1.D0
          DRELPR = D1MACH(4)
          CHECKL = .FALSE.
          FILTER = .TRUE.
          LENX = 2* (NCOLS+MCON) + 2
          ISCALE = 1
          IGO = 1
          ACCUM = .FALSE.
          PRETRI = .TRUE.
          LOPT = 0
          LP = 0
          LDS = 0
C     DO FOREVER
   30     CONTINUE
          LP = LP + LDS
          IP = IOPT(LP+1)
          JP = ABS(IP)
C
C     TEST FOR NO MORE OPTIONS TO CHANGE.
          IF (IP.EQ.99) THEN
              IF (LOPT.EQ.0) LOPT = LP+1
C
C     SEND COL. SCALING TO DBOLS().
              IDOPE(4)=1
C
C     NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
              IDOPE(5)=1
C
C     CHANGE PRETRIANGULARIZATION FACTOR IN DBOLSM().
              IDOPE(1) = NCOLS + MCON + 1
C
C     PASS WEIGHT TO DBOLSM() FOR RANK TEST.
              IDOPE(2) = NCOLS + MCON + 2
              IDOPE(3) = MCON
C     EXIT FOREVER
              GO TO 50
          ELSE IF (JP.EQ.99) THEN
              LDS = 1
C     CYCLE FOREVER
              GO TO 50
          ELSE IF (JP.EQ.1) THEN
              IF (IP.GT.0) THEN
C
C     SET UP DIRECTION FLAG LOCATION, ROW STACKING POINTER
C     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
                  LOCACC = LP + 2
C
C                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
C     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
C                  IOPT(LOCACC+1)=ROW STACKING POINTER.
C                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
C     USER ACTION WITH THIS OPTION..
C      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).)
C      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
C       ROW OF W(*,*) BELOW THE ROWS FOR THE CONSTRAINT MATRIX C.
C       SET IOPT(LOCACC+2)=NO. OF LEAST SQUARES EQUATIONS IN BLOCK.
C              LOOP
C              CALL DBOCLS()
C
C                  IF(IOPT(LOCACC) .EQ. 1) THEN
C                      STACK EQUAS. INTO W(*,*), STARTING AT
C                      ROW IOPT(LOCACC+1).
C                       INTO W(*,*).
C                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
C                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
C                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN
C                      (PROCESS IS OVER. EXIT LOOP.)
C                  ELSE
C                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
C                  END IF
C              END LOOP
                  IOPT(LOCACC+1) = MCON + 1
                  ACCUM = .TRUE.
                  IOPT(LOCACC) = IGO
              END IF
              LDS = 4
C     CYCLE FOREVER
              GO TO 30
          ELSE IF (JP.EQ.2) THEN
              IF (IP.GT.0) THEN
C
C     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
                  LOCDIM = LP + 2
C
C     LMDW.GE.MCON+MAX(MOUT,NCOLS), IF MCON.GT.0 .AND FILTER
C     LMDW.GE.MCON+MOUT, OTHERWISE
C
C     LNDW.GE.NCOLS+MCON+1
C     LLB .GE.NCOLS+MCON
C     LLX .GE.2*(NCOLS+MCON)+2+EXTRA REQD. IN OPTIONS.
C     LLRW.GE.6*NCOLS+5*MCON
C     LLIW.GE.2*(NCOLS+MCON)
C     LIOP.GE. AMOUNT REQD. FOR OPTION ARRAY.
                  LMDW = IOPT(LOCDIM)
                  LNDW = IOPT(LOCDIM+1)
                  LLB = IOPT(LOCDIM+2)
                  LLX = IOPT(LOCDIM+3)
                  LLRW = IOPT(LOCDIM+4)
                  LLIW = IOPT(LOCDIM+5)
                  LIOPT = IOPT(LOCDIM+6)
                  CHECKL = .TRUE.
              END IF
              LDS = 8
C     CYCLE FOREVER
              GO TO 30
C
C     OPTION TO MODIFY THE COLUMN SCALING.
          ELSE IF (JP.EQ.3) THEN
              IF (IP.GT.0) THEN
                  ISCALE = IOPT(LP+2)
C
C     SEE THAT ISCALE IS 1 THRU 3.
                  IF (ISCALE.LT.1 .OR. ISCALE.GT.3) THEN
                      NERR = 48
                      NCHAR = 41
                      CALL XERRWV(
     *                       'DBOCLS(). ISCALE OPTION=(I1) MUST BE 1-3.'
     *                            ,NCHAR,NERR,LEVEL,1,ISCALE,IDUM,0,
     *                            RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 260
                  END IF
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
C     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
          ELSE IF (JP.EQ.4) THEN
              IF (IP.GT.0) THEN
                  ISCALE = 4
                  IF (IOPT(LP+2).LE.0) THEN
                      NERR = 49
                      NCHAR = 86
                      CALL XERRWV(
     *'DBOCLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED COLUMN SCAL
     *ING MUST BE POSITIVE.',NCHAR,NERR,LEVEL,1,IOPT(LP+2),IDUM,0,RDUM,
     *                            RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 260
                  END IF
                  CALL DCOPY(NCOLS,X(NCOLS+IOPT(LP+2)),1,RW,1)
                  LENX = LENX + NCOLS
                  DO 40 J = 1,NCOLS
                      IF (RW(J).LE.ZERO) THEN
                          NERR = 50
                          NCHAR = 84
                          CALL XERRWV(
     *'DBOCLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE. COMP.
     * (I1)   NOW = (R1).',NCHAR,NERR,LEVEL,1,J,IDUM,1,RW(J),RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                          GO TO 260
                      END IF
   40             CONTINUE
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLS().
          ELSE IF (JP.EQ.5) THEN
              IF (IP.GT.0) THEN
                  LOPT = IOPT(LP+2)
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
C     (NO LONGER USED.) OPTION NOW MUST BE PASSED IMBEDDED IN
C     OPTION ARRAY FOR DBOLS().
          ELSE IF (JP.EQ.6) THEN
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS A
C     POINTER VALUE TO SKIP TO NEXT.
          ELSE IF (JP.EQ.7) THEN
              IF (IP.GT.0) THEN
                  LP = IOPT(LP+2)-1
                  LDS = 0
              ELSE
                  LDS = 2
              END IF
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION AVOIDS THE CONSTRAINT RESOLVING PHASE FOR
C     THE LINEAR CONSTRAINTS C*X=Y.
          ELSE IF (JP.EQ.8) THEN
              FILTER = .NOT. (IP.GT.0)
              LDS = 1
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION SUPPRESSES PRETRIANGULARIZATION OF THE LEAST
C     SQUARES EQATIONS.
          ELSE IF (JP.EQ.9) THEN
              PRETRI = .NOT. (IP.GT.0)
              LDS = 1
C     CYCLE FOREVER
              GO TO 30
C
C     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
          ELSE
              NERR = 51
              NCHAR = 48
              CALL XERRWV(
     *                'DBOCLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED.'
     *                    ,NCHAR,NERR,LEVEL,1,JP,IDUM,0,IDUM,IDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 260
          END IF
C     END FOREVER
C     END PROCEDURE
   50     CONTINUE
          IF (CHECKL) THEN
C     DO(CHECK LENGTHS OF ARRAYS)
C     PROCEDURE(CHECK LENGTHS OF ARRAYS)
C
C     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
C     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
           IF(FILTER .AND. .NOT.ACCUM) THEN
                MDWL=MCON+MAX(MROWS,NCOLS)
           ELSE IF (ACCUM) THEN
                MDWL=MCON+NCOLS+1
           ELSE
                MDWL=MCON+NCOLS
           END IF
              IF (LMDW.LT.MDWL) THEN
                  NERR = 41
                  NCHAR = 88
                  CALL XERRWV(
     *'DBOCLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE NUMBER
     *OF EFFECTIVE ROWS=(I2).',NCHAR,NERR,LEVEL,2,LMDW,
     *MDWL,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LNDW.LT.NCOLS+MCON+1) THEN
                  NERR = 42
                  NCHAR = 75
                  CALL XERRWV(
     *'DBOCLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+MC
     *ON+1=(I2).',NCHAR,NERR,LEVEL,2,LNDW,NCOLS+MCON+1,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LLB.LT.NCOLS+MCON) THEN
                  NERR = 43
                  NCHAR = 94
                  CALL XERRWV(
     *'DBOCLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1)
     *MUST BE .GE. NCOLS+MCON=(I2).',NCHAR,NERR,LEVEL,2,LLB,NCOLS+MCON,
     *                        0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LLX.LT.LENX) THEN
                  NERR = 44
                  NCHAR = 71
                  CALL XERRWV(
     *'DBOCLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. THE REQD. LENGTH
     *=(I2).',NCHAR,NERR,LEVEL,2,LLX,LENX,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LLRW.LT.6*NCOLS+5*MCON) THEN
                  NERR = 45
                  NCHAR = 70
                  CALL XERRWV(
     *'DBOCLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 6*NCOLS+5*MCON=
     *(I2).',NCHAR,NERR,LEVEL,2,LLRW,6*NCOLS+5*MCON,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LLIW.LT.2*NCOLS+2*MCON) THEN
                  NERR = 46
                  NCHAR = 69
                  CALL XERRWV(
     *'DBOCLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS+2*MCON=(
     *I2).',NCHAR,NERR,LEVEL,2,LLIW,2*NCOLS+2*MCON,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
              IF (LIOPT.LT.LP+17) THEN
                  NERR = 47
                  NCHAR = 72
                  CALL XERRWV(
     *'DBOCLS(). THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. LEN
     *.=(I2).',NCHAR,NERR,LEVEL,2,LIOPT,LP+17,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 260
              END IF
C     END PROCEDURE
          END IF
      END IF
C
C     OPTIONALLY GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
C     EQUATIONS AND DIRECTIONS FOR PROCESSING THESE EQUATIONS.
C     DO(ACCUMULATE LEAST SQUARES EQUATIONS)
C     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS)
      IF (ACCUM) THEN
          MROWS = IOPT(LOCACC+1) - 1 - MCON
          INROWS = IOPT(LOCACC+2)
          MNEW = MROWS + INROWS
          IF (MNEW.LT.0 .OR. MNEW+MCON.GT.MDW) THEN
              NERR = 52
              NCHAR = 66
              CALL XERRWV(
     *'DBOCLS(). NO. OF ROWS=(I1) MUST BE .GE. 0 .AND. .LE.MDW-MCON=(I2)
     *',NCHAR,NERR,LEVEL,2,MNEW,MDW-MCON,0,RDUM,RDUM)
C    (RETURN TO USER PROGRAM UNIT)
              GO TO 260
          END IF
      END IF
C
C     USE THE SOFTWARE OF DBOLS( ) FOR THE TRIANGULARIZATION OF THE
C     LEAST SQUARES MATRIX.  THIS MAY INVOLVE A SYSTALTIC INTERCHANGE
C     OF PROCESSING POINTERS BETWEEN THE CALLING AND CALLED (DBOLS())
C     PROGRAM UNITS.
      JOPT(01) = 1
      JOPT(02) = 2
      JOPT(04) = MROWS
      JOPT(05) = 99
      IRW = NCOLS + 1
      IIW = 1
      IF (ACCUM .OR. PRETRI) THEN
C
C     NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
              IDOPE(5)=0
          CALL DBOLS(W(MCON+1,1),MDW,MOUT,NCOLS,BL,BU,IND,JOPT,X,RNORM,
     *               MODE,RW(IRW),IW(IIW))
      ELSE
          MOUT = MROWS
      END IF
      IF (ACCUM) THEN
          ACCUM = IOPT(LOCACC) .EQ. 1
          IOPT(LOCACC+1) = JOPT(03) + MCON
          MROWS = MIN(NCOLS+1,MNEW)
      END IF
C     END PROCEDURE
      IF (ACCUM) RETURN
C     DO(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM)
C     PROCEDURE(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM)
C
C     MOVE RIGHT HAND SIDE OF LEAST SQUARES EQUATIONS.
      CALL DCOPY(MOUT,W(MCON+1,NCOLS+1),1,W(MCON+1,NCOLS+MCON+1),1)
      IF (MCON.GT.0 .AND. FILTER) THEN
C
C     PROJECT THE LINEAR CONSTRAINTS INTO A REACHABLE SET.
          DO 60 I = 1,MCON
              CALL DCOPY(NCOLS,W(I,1),MDW,W(MCON+1,NCOLS+I),1)
   60     CONTINUE
C
C      PLACE (-)IDENTITY MATRIX AFTER CONSTRAINT DATA.
          DO 70 J = NCOLS + 1,NCOLS + MCON + 1
              W(1,J) = ZERO
              CALL DCOPY(MCON,W(1,J),0,W(1,J),1)
   70     CONTINUE
          W(1,NCOLS+1) = -ONE
          CALL DCOPY(MCON,W(1,NCOLS+1),0,W(1,NCOLS+1),MDW+1)
C
C     OBTAIN A 'FEASIBLE POINT' FOR THE LINEAR CONSTRAINTS.
          JOPT(01) = 99
          IRW = NCOLS + 1
          IIW = 1
C
C     NOTE THAT DBOLS() WAS CALLED BY DBOCLS()
              IDOPE(5)=0
          CALL DBOLS(W,MDW,MCON,NCOLS+MCON,BL,BU,IND,JOPT,X,RNORMC,
     *               MODEC,RW(IRW),IW(IIW))
C
C     ENLARGE THE BOUNDS SET, IF REQUIRED, TO INCLUDE POINTS THAT
C     CAN BE REACHED.
          DO 130 J = NCOLS + 1,NCOLS + MCON
              ICASE = IND(J)
              IF (ICASE.LT.4) THEN
                  T = DDOT(NCOLS,W(MCON+1,J),1,X,1)
              END IF
              GO TO (80,90,100,110),ICASE
              GO TO 120
C     CASE 1
   80         BL(J) = MIN(T,BL(J))
              GO TO 120
C     CASE 2
   90         BU(J) = MAX(T,BU(J))
              GO TO 120
C     CASE 3
  100         BL(J) = MIN(T,BL(J))
              BU(J) = MAX(T,BU(J))
              GO TO 120
C     CASE 4
  110         CONTINUE
  120         CONTINUE
  130     CONTINUE
C
C     MOVE CONSTRAINT DATA BACK TO THE ORIGINAL AREA.
          DO 140 J = NCOLS + 1,NCOLS + MCON
              CALL DCOPY(NCOLS,W(MCON+1,J),1,W(J-NCOLS,1),MDW)
  140     CONTINUE
      END IF
      IF (MCON.GT.0) THEN
          DO 150 J = NCOLS + 1,NCOLS + MCON
              W(MCON+1,J) = ZERO
              CALL DCOPY(MOUT,W(MCON+1,J),0,W(MCON+1,J),1)
  150     CONTINUE
C
C     PUT IN (-)IDENTITY MATRIX (POSSIBLY) ONCE AGAIN.
          DO 160 J = NCOLS + 1,NCOLS + MCON + 1
              W(1,J) = ZERO
              CALL DCOPY(MCON,W(1,J),0,W(1,J),1)
  160     CONTINUE
          W(1,NCOLS+1) = -ONE
          CALL DCOPY(MCON,W(1,NCOLS+1),0,W(1,NCOLS+1),MDW+1)
      END IF
C
C     COMPUTE NOMINAL COLUMN SCALING FOR THE UNWEIGHTED MATRIX.
      CNORM = ZERO
      ANORM = ZERO
      DO 170 J = 1,NCOLS
          T1 = DASUM(MCON,W(1,J),1)
          T2 = DASUM(MOUT,W(MCON+1,1),1)
          T = T1 + T2
          IF (T.EQ.ZERO) T = ONE
          CNORM = MAX(CNORM,T1)
          ANORM = MAX(ANORM,T2)
          X(NCOLS+MCON+J) = ONE/T
  170 CONTINUE
      GO TO (180,190,210,220),ISCALE
      GO TO 230
C     CASE 1
  180 CONTINUE
      GO TO 230
C     CASE 2
C
C     SCALE COLS. (BEFORE WEIGHTING) TO HAVE LENGTH ONE.
  190 DO 200 J = 1,NCOLS
          T = DNRM2(MCON+MOUT,W(1,J),1)
          IF (T.EQ.ZERO) T = ONE
          X(NCOLS+MCON+J) = ONE/T
  200 CONTINUE
      GO TO 230
C     CASE 3
C
C     SUPPRESS SCALING (USE UNIT MATRIX).
  210 X(NCOLS+MCON+1) = ONE
      CALL DCOPY(NCOLS,X(NCOLS+MCON+1),0,X(NCOLS+MCON+1),1)
      GO TO 230
C     CASE 4
C
C     THE USER HAS PROVIDED SCALING.
  220 CALL DCOPY(NCOLS,RW,1,X(NCOLS+MCON+1),1)
  230 CONTINUE
      DO 240 J = NCOLS + 1,NCOLS + MCON
          X(NCOLS+MCON+J) = ONE
  240 CONTINUE
C
C     WEIGHT THE LEAST SQUARES EQUATIONS.
      WT = SQRT(DRELPR)
      IF (ANORM.GT.ZERO) WT = WT/ANORM
      IF (CNORM.GT.ZERO) WT = WT*CNORM
      DO 250 I = 1,MOUT
          CALL DSCAL(NCOLS,WT,W(I+MCON,1),MDW)
  250 CONTINUE
      CALL DSCAL(MOUT,WT,W(MCON+1,MCON+NCOLS+1),1)
      LRW = 1
      LIW = 1
C
C     SET THE NEW TRIANGULARIZATION FACTOR.
      X(NCOLS+MCON+IDOPE(1))= ZERO
C
C     SET THE WEIGHT TO USE IN COMPONENTS .GT. MCON,
C     WHEN MAKING LINEAR INDEPENDENCE TEST.
      X(NCOLS+MCON+IDOPE(2))= ONE/WT
      IDOPE(5)=1
      CALL DBOLS(W,MDW,MOUT+MCON,NCOLS+MCON,BL,BU,IND,IOPT(LOPT),X,
     *           RNORM,MODE,RW(LRW),IW(LIW))
C     END PROCEDURE
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  260 IF(MODE.GE.0)MODE=-NERR
      IGO = 0
      RETURN
C     END PROGRAM
      END
      SUBROUTINE DBOLS(W,MDW,MROWS,NCOLS,BL,BU,IND,IOPT,X,RNORM,MODE,
     *   RW,IW)
C***BEGIN PROLOGUE  DBOLS
C***DATE WRITTEN   821220   (YYMMDD)
C***REVISION DATE  880722   (yymmdd)
C***CATEGORY NO.  K1A2A,G2E,G2H1,G2H2
C***KEYWORDS  BOUNDS,CONSTRAINTS,INEQUALITY,LEAST SQUARES,LINEAR
C***AUTHOR  HANSON, R. J., SNLA
C***PURPOSE  Solve the problem
C                 E*X = F (in the least  squares  sense)
C            with bounds on selected X values.
C***DESCRIPTION
C
C     The user must have dimension statements of the form:
C
C       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
C      * X(NCOLS+NX), RW(5*NCOLS)
C       INTEGER IND(NCOLS), IOPT(1+NI), IW(2*NCOLS)
C
C     (here NX=number of extra locations required for option 4; NX=0
C     for no options; NX=NCOLS if this option is in use. Here NI=number
C     of extra locations required for options 1-6; NI=0 for no
C     options.)
C
C   INPUT
C   -----
C
C    --------------------
C    W(MDW,*),MROWS,NCOLS
C    --------------------
C     The array W(*,*) contains the matrix [E:F] on entry. The matrix
C     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
C     the array W(*,*) with E occupying the first NCOLS columns and the
C     right side vector F in column NCOLS+1. The row dimension, MDW, of
C     the array W(*,*) must satisfy the inequality MDW .ge. MROWS.
C     Other values of MDW are errrors. The values of MROWS and NCOLS
C     must be positive. Other values are errors. There is an exception
C     to this when using option 1 for accumulation of blocks of
C     equations. In that case MROWS is an OUTPUT variable ONLY, and the
C     matrix data for [E:F] is placed in W(*,*), one block of rows at a
C     time.  MROWS contains the number of rows in the matrix after
C     triangularizing several blocks of equations. This is an OUTPUT
C     parameter ONLY when option 1 is used. See IOPT(*) CONTENTS
C     for details about option 1.
C
C    ------------------
C    BL(*),BU(*),IND(*)
C    ------------------
C     These arrays contain the information about the bounds that the
C     solution values are to satisfy. The value of IND(J) tells the
C     type of bound and BL(J) and BU(J) give the explicit values for
C     the respective upper and lower bounds.
C
C    1.    For IND(J)=1, require X(J) .ge. BL(J).
C          (the value of BU(J) is not used.)
C    2.    For IND(J)=2, require X(J) .le. BU(J).
C          (the value of BL(J) is not used.)
C    3.    For IND(J)=3, require X(J) .ge. BL(J) and
C                                X(J) .le. BU(J).
C    4.    For IND(J)=4, no bounds on X(J) are required.
C          (the values of BL(J) and BU(J) are not used.)
C
C     Values other than 1,2,3 or 4 for IND(J) are errors. In the case
C     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J)
C     is an error.
C
C    -------
C    IOPT(*)
C    -------
C     This is the array where the user can specify nonstandard options
C     for DBOLSM( ). Most of the time this feature can be ignored by
C     setting the input value IOPT(1)=99. Occasionally users may have
C     needs that require use of the following subprogram options. For
C     details about how to use the options see below: IOPT(*) CONTENTS.
C
C     Option Number   Brief Statement of Purpose
C     ------ ------   ----- --------- -- -------
C           1         Return to user for accumulation of blocks
C                     of least squares equations.
C           2         Check lengths of all arrays used in the
C                     subprogram.
C           3         Standard scaling of the data matrix, E.
C           4         User provides column scaling for matrix E.
C           5         Provide option array to the low-level
C                     subprogram DBOLSM( ).
C           6         Move the IOPT(*) processing pointer.
C           7         User has called DBOLS() directly.
C          99         No more options to change.
C
C    ----
C    X(*)
C    ----
C     This array is used to pass data associated with option 4. Ignore
C     this parameter if this option is not used. Otherwise see below:
C     IOPT(*) CONTENTS.
C
C    OUTPUT
C    ------
C
C    ----------
C    X(*),RNORM
C    ----------
C     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
C     the constrained least squares problem. The value RNORM is the
C     minimum residual vector length.
C
C    ----
C    MODE
C    ----
C     The sign of MODE determines whether the subprogram has completed
C     normally, or encountered an error condition or abnormal status. A
C     value of MODE .ge. 0 signifies that the subprogram has completed
C     normally. The value of MODE (.GE. 0) is the number of variables
C     in an active status: not at a bound nor at the value ZERO, for
C     the case of free variables. A negative value of MODE will be one
C     of the cases -37,-36,...,-22, or -17,...,-2. Values .lt. -1
C     correspond to an abnormal completion of the subprogram. To
C     understand the abnormal completion codes see below: ERROR
C     MESSAGES for DBOLS( ). AN approximate solution will be returned
C     to the user only when max. iterations is reached, MODE=-22.
C     Values for MODE=-37,...,-22 come from the low-level subprogram
C     DBOLSM(). See the section ERROR MESSAGES for DBOLSM() in the
C     documentation for DBOLSM().
C
C    -----------
C    RW(*),IW(*)
C    -----------
C     These are working arrays with 5*NCOLS and 2*NCOLS entries.
C     (normally the user can ignore the contents of these arrays,
C     but they must be dimensioned properly.)
C
C    IOPT(*) CONTENTS
C    ------- --------
C     The option array allows a user to modify internal variables in
C     the subprogram without recompiling the source code. A central
C     goal of the initial software design was to do a good job for most
C     people. Thus the use of options will be restricted to a select
C     group of users. The processing of the option array proceeds as
C     follows: a pointer, here called LP, is initially set to the value
C     1. This value is updated as each option is processed. At the
C     pointer position the option number is extracted and used for
C     locating other information that allows for options to be changed.
C     The portion of the array IOPT(*) that is used for each option is
C     fixed; the user and the subprogram both know how many locations
C     are needed for each option. A great deal of error checking is
C     done by the subprogram on the contents of the option array.
C     Nevertheless it is still possible to give the subprogram optional
C     input that is meaningless. For example option 4 uses the
C     locations X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing
C     scaling data. The user must manage the allocation of these
C     locations.
C
C   1
C   -
C     This option allows the user to solve problems with a large number
C     of rows compared to the number of variables. The idea is that the
C     subprogram returns to the user (perhaps many times) and receives
C     new least squares equations from the calling program unit.
C     Eventually the user signals "that's all" and then computes the
C     solution with one final call to subprogram DBOLS( ). The value of
C     MROWS is an OUTPUT variable when this option is used. Its value
C     is always in the range 0 .le. MROWS .le. NCOLS+1. It is equal to
C     the number of rows after the triangularization of the entire set
C     of equations. If LP is the processing pointer for IOPT(*), the
C     usage for the sequential processing of blocks of equations is
C
C        IOPT(LP)=1
C        Move block of equations to W(*,*) starting at
C        the first row of W(*,*).
C        IOPT(LP+3)=# of rows in the block; user defined
C
C     The user now calls DBOLS( ) in a loop. The value of IOPT(LP+1)
C     directs the user's action. The value of IOPT(LP+2) points to
C     where the subsequent rows are to be placed in W(*,*).
C
C      .<LOOP
C      . CALL DBOLS()
C      . IF(IOPT(LP+1) .EQ. 1) THEN
C      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED
C      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN
C      .    W(*,*) STARTING AT ROW IOPT(LP+2).
C      .
C      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN
C      .       IOPT(LP+1)=2
C      .<------CYCLE LOOP
C      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN
C      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0
C      . ELSE
C      . ERROR CONDITION; SHOULD NOT HAPPEN.
C      .<END LOOP
C
C     Use of this option adds 4 to the required length of IOPT(*).
C
C
C   2
C   -
C     This option is useful for checking the lengths of all arrays used
C     by DBOLS() against their actual requirements for this problem.
C     The idea is simple: the user's program unit passes the declared
C     dimension information of the arrays. These values are compared
C     against the problem-dependent needs within the subprogram. If any
C     of the dimensions are too small an error message is printed and a
C     negative value of MODE is returned, -11 to -17. The printed error
C     message tells how long the dimension should be. If LP is the
C     processing pointer for IOPT(*),
C
C        IOPT(LP)=2
C        IOPT(LP+1)=Row dimension of W(*,*)
C        IOPT(LP+2)=Col. dimension of W(*,*)
C        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*)
C        IOPT(LP+4)=Dimension of X(*)
C        IOPT(LP+5)=Dimension of RW(*)
C        IOPT(LP+6)=Dimension of IW(*)
C        IOPT(LP+7)=Dimension of IOPT(*)
C         .
C        CALL DBOLS()
C
C     Use of this option adds 8 to the required length of IOPT(*).
C
C   3
C   -
C     This option changes the type of scaling for the data matrix E.
C     Nominally each nonzero column of E is scaled so that the
C     magnitude of its largest entry is equal to the value ONE. If LP
C     is the processing pointer for IOPT(*),
C
C        IOPT(LP)=3
C        IOPT(LP+1)=1,2 or 3
C            1= Nominal scaling as noted;
C            2= Each nonzero column scaled to have length ONE;
C            3= Identity scaling; scaling effectively suppressed.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   4
C   -
C     This option allows the user to provide arbitrary (positive)
C     column scaling for the matrix E. If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=4
C        IOPT(LP+1)=IOFF
C        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1)
C        = Positive scale factors for cols. of E.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*) and
C     NCOLS to the required length of X(*).
C
C   5
C   -
C     This option allows the user to provide an option array to the
C     low-level subprogram DBOLSM(). If LP is the processing pointer
C     for IOPT(*),
C
C        IOPT(LP)=5
C        IOPT(LP+1)= Position in IOPT(*) where option array
C                    data for DBOLSM() begins.
C         .
C        CALL DBOLS()
C
C     Use of this option adds 2 to the required length of IOPT(*).
C
C   6
C   -
C     Move the processing pointer (either forward or backward) to the
C     location IOPT(LP+1). The processing point is moved to entry
C     LP+2 of IOPT(*) if the option is left with -6 in IOPT(LP).  For
C     example to skip over locations 3,...,NCOLS+2 of IOPT(*),
C
C       IOPT(1)=6
C       IOPT(2)=NCOLS+3
C       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
C       IOPT(NCOLS+3)=99
C       CALL DBOLS()
C
C     CAUTION: Misuse of this option can yield some very hard
C     -to-find bugs.  Use it with care.
C
C   7
C   -
C     If the user is calling DBOLS() directly, use this option.
C     (This is necessary because DBOCLS() uses DBOLS() as a
C     low-level subprogram.  Due to weighting required within
C     DBOCLS(), the two cases must be known.) For example,
C
C       IOPT(1)=7
C       IOPT(1)=99
C
C   99
C   --
C     There are no more options to change.
C
C     Only option numbers -99, -7,-6,-5,...,-1, 1,2,...,7, and 99 are
C     permitted. Other values are errors. Options -99,-1,...,-7 mean
C     that the repective options 99,1,...,7 are left at their default
C     values. An example is the option to modify the (rank) tolerance:
C
C       IOPT(1)=-3 Option is recognized but not changed
C       IOPT(2)=2  Scale nonzero cols. to have length ONE
C       IOPT(3)=99
C
C    ERROR MESSAGES for DBOLS()
C    ----- -------- --- -------
C
C WARNING IN...
C DBOLS(). MDW=(I1) MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         2
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         3
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         0
C ERROR NUMBER =         4
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, R1=    0.
C           IN ABOVE MESSAGE, R2=    ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         6
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). ISCALE OPTION=(I1) MUST BE 1-3.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         7
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED  COLUMN SCALING
C MUST BE POSITIVE.
C           IN ABOVE MESSAGE, I1=         0
C ERROR NUMBER =         8
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE.
C COMPONENT (I1) NOW = (R1).
C           IN ABOVE MESSAGE, I1=        ND. .LE. MDW=(I2).
C           IN ABOVE MESSAGE, I1=         1
C           IN ABOVE MESSAGE, I2=         0
C ERROR NUMBER =        10
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS().THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE.THE NUMBER OF ROWS=
C (I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        11
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+1=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        12
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS().THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) MUST BE
C .GE. NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        13
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. THE REQD. LENGTH=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        14
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 5*NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         3
C ERROR NUMBER =        15
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         2
C ERROR NUMBER =        16
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C WARNING IN...
C DBOLS() THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. LEN.=(I2).
C           IN ABOVE MESSAGE, I1=         0
C           IN ABOVE MESSAGE, I2=         1
C ERROR NUMBER =        17
C (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C***REFERENCES  HANSON, R. J. LINEAR LEAST SQUARES WITH BOUNDS AND
C                 LINEAR CONSTRAINTS, SIAM J. SCI. STAT. COMPUT., VOL. 7,
C                 NO. 3, JULY, 1986.
C***ROUTINES CALLED  IDAMAX,DBOLSM,DCOPY,DNRM2,DROT,DROTG,XERRWV
C***COMMON BLOCKS    DBOCOM
C***END PROLOGUE  DBOLS
C
C     SOLVE LINEAR LEAST SQUARES SYSTEM WITH BOUNDS ON
C     SELECTED VARIABLES.
C     REVISED 880722-1100
C     REVISED YYMMDD-HHMM
C     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
C     EDITING AT THE CARD 'C++'.
C     CHANGE THIS SUBPROGRAM NAME TO DBOLS AND THE STRINGS
C     /SCOPY/ TO /DCOPY/, /DBOL/ TO /DBOL/,
C     /SNRM2/ TO /DNRM2/, /ISAMAX/ TO /IDAMAX/,
C     /SROTG/ TO /DROTG/, /SROT/ TO /DROT/, /E0/ TO /D0/,
C     /REAL            / TO /DOUBLE PRECISION/.
C ++
      DOUBLE PRECISION W(MDW,*),BL(*),BU(*),X(*),RW(*)
      DOUBLE PRECISION SC, SS, ONE, DNRM2, RNORM, ZERO
C
C     THIS VARIABLE SHOULD REMAIN TYPE REAL.
      REAL RDUM
      INTEGER IND(*),IOPT(*),IW(*)
      LOGICAL CHECKL
      COMMON /DBOCOM/ IDOPE(5)
      SAVE IGO,LOCACC,LOPT,ISCALE
      DATA IGO/0/
C***FIRST EXECUTABLE STATEMENT  DBOLS
      NERR = 0
      MODE = 0
      LEVEL = 1
      IF (IGO.EQ.0) THEN
C     DO(CHECK VALIDITY OF INPUT DATA)
C     PROCEDURE(CHECK VALIDITY OF INPUT DATA)
C
C     SEE THAT MDW IS .GT.0. GROSS CHECK ONLY.
          IF (MDW.LE.0) THEN
              NERR = 2
              NCHAR = 35
              CALL XERRWV('DBOLS(). MDW=(I1) MUST BE POSITIVE.',NCHAR,
     .                    NERR,LEVEL,1,MDW,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          END IF
C
C     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE.
          IF (NCOLS.LE.0) THEN
              NERR = 3
              NCHAR = 58
              CALL XERRWV(
     .      'DBOLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE.'
     .                    ,NCHAR,NERR,LEVEL,1,NCOLS,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          END IF
C
C     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED.
          DO 10 J = 1,NCOLS
              IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
                  NERR = 4
                  NCHAR = 45
                  CALL XERRWV(
     .                   'DBOLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.'
     .                        ,NCHAR,NERR,LEVEL,2,J,IND(J),0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
   10     CONTINUE
C
C     SEE THAT BOUNDS ARE CONSISTENT.
          DO 20 J = 1,NCOLS
              IF (IND(J).EQ.3) THEN
                  IF (BL(J).GT.BU(J)) THEN
                      NERR = 5
                      NCHAR = 57
                      CALL XERRWV(
     .       'DBOLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2).'
     .                            ,NCHAR,NERR,LEVEL,1,J,IDUM,2,BL(J),
     .                            BU(J))
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  END IF
              END IF
   20     CONTINUE
C     END PROCEDURE
C     DO(PROCESS OPTION ARRAY)
C     PROCEDURE(PROCESS OPTION ARRAY)
          ZERO = 0.D0
          ONE = 1.D0
          CHECKL = .FALSE.
          LENX = NCOLS
          ISCALE = IDOPE(4)
          IGO = 2
          LOPT = 0
          LP = 0
          LDS = 0
   30     CONTINUE
          LP = LP + LDS
          IP = IOPT(LP+1)
          JP = ABS(IP)
C
C     TEST FOR NO MORE OPTIONS.
          IF (IP.EQ.99) THEN
              IF (LOPT.EQ.0) LOPT = LP + 1
              GO TO 50
          ELSE IF (JP.EQ.99) THEN
              LDS = 1
              GO TO 30
          ELSE IF (JP.EQ.1) THEN
              IF (IP.GT.0) THEN
C
C     SET UP DIRECTION FLAG, ROW STACKING POINTER
C     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS.
                  LOCACC = LP + 2
C
C                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION.
C     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2.
C                  IOPT(LOCACC+1)=ROW STACKING POINTER.
C                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS.
C     USER ACTION WITH THIS OPTION..
C      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).
C      MUST ALSO START PROCESS WITH IOPT(LOCACC)=1.)
C      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST
C       ROW OF W(*,*).  SET IOPT(LOCACC+2)=NO. OF ROWS IN BLOCK.)
C              LOOP
C              CALL DBOLS()
C
C                  IF(IOPT(LOCACC) .EQ. 1) THEN
C                      STACK EQUAS., STARTING AT ROW IOPT(LOCACC+1),
C                       INTO W(*,*).
C                       SET IOPT(LOCACC+2)=NO. OF EQUAS.
C                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2.
C                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN
C                      (PROCESS IS OVER. EXIT LOOP.)
C                  ELSE
C                      (ERROR CONDITION. SHOULD NOT HAPPEN.)
C                  END IF
C              END LOOP
C              SET IOPT(LOCACC-1)=-OPTION NUMBER FOR SEQ. ACCUMULATION.
C              CALL DBOLS( )
                  IOPT(LOCACC+1) = 1
                  IGO = 1
              END IF
              LDS = 4
              GO TO 30
          ELSE IF (JP.EQ.2) THEN
              IF (IP.GT.0) THEN
C
C     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS.
                  LOCDIM = LP + 2
C
C     LMDW.GE.MROWS
C     LNDW.GE.NCOLS+1
C     LLB .GE.NCOLS
C     LLX .GE.NCOLS+EXTRA REQD. IN OPTIONS.
C     LLRW.GE.5*NCOLS
C     LLIW.GE.2*NCOLS
C     LIOP.GE. AMOUNT REQD. FOR IOPTION ARRAY.
                  LMDW = IOPT(LOCDIM)
                  LNDW = IOPT(LOCDIM+1)
                  LLB = IOPT(LOCDIM+2)
                  LLX = IOPT(LOCDIM+3)
                  LLRW = IOPT(LOCDIM+4)
                  LLIW = IOPT(LOCDIM+5)
                  LIOPT = IOPT(LOCDIM+6)
                  CHECKL = .TRUE.
              END IF
              LDS = 8
              GO TO 30
C
C     OPTION TO MODIFY THE COLUMN SCALING.
          ELSE IF (JP.EQ.3) THEN
              IF (IP.GT.0) THEN
                  ISCALE = IOPT(LP+2)
C
C     SEE THAT ISCALE IS 1 THRU 3.
                  IF (ISCALE.LT.1 .OR. ISCALE.GT.3) THEN
                      NERR = 7
                      NCHAR = 40
                      CALL XERRWV(
     .                        'DBOLS(). ISCALE OPTION=(I1) MUST BE 1-3.'
     .                            ,NCHAR,NERR,LEVEL,1,ISCALE,IDUM,0,
     .                            RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  END IF
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE
C     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)).
          ELSE IF (JP.EQ.4) THEN
              IF (IP.GT.0) THEN
                  ISCALE = 4
                  IF (IOPT(LP+2).LE.0) THEN
                      NERR = 8
                      NCHAR = 85
                      CALL XERRWV(
     .'DBOLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED COLUMN SCALI
     .NG MUST BE POSITIVE.',NCHAR,NERR,LEVEL,1,IOPT(LP+2),IDUM,0,RDUM,
     .                            RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                      GO TO 190
                  END IF
                  CALL DCOPY(NCOLS,X(NCOLS+IOPT(LP+2)),1,RW,1)
                  LENX = LENX + NCOLS
                  DO 40 J = 1,NCOLS
                      IF (RW(J).LE.ZERO) THEN
                          NERR = 9
                          NCHAR = 85
                          CALL XERRWV(
     .'DBOLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE. COMPON
     .ENT (I1) NOW = (R1).',NCHAR,NERR,LEVEL,1,J,IDUM,1,RW(J),RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                          GO TO 190
                      END IF
   40             CONTINUE
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM().
          ELSE IF (JP.EQ.5) THEN
              IF (IP.GT.0) THEN
                  LOPT = IOPT(LP+2)
              END IF
              LDS = 2
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS THE PLACE TO
C     MOVE AT THE NEXT STEP OF PROCESSESING.
          ELSE IF (JP.EQ.6) THEN
              IF (IP.GT.0) THEN
                  LP = IOPT(LP+2)-1
                  LDS = 0
              ELSE
                  LDS = 2
              END IF
C     CYCLE FOREVER
              GO TO 30
C
C     THIS OPTION PROVIDES INFORMATION ABOUT WHO CALLED DBOLS.
C     IT WAS EITHER DBOCLS() OR THE USER.
          ELSE IF (JP.EQ.7) THEN
              LDS=1
              IF (IP.GT.0) THEN
                IDOPE(5)=0
                ISCALE=1
              END IF
C     CYCLE FOREVER
              GO TO 30
C
C     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION.
          ELSE
              NERR = 6
              NCHAR = 47
              CALL XERRWV(
     .                 'DBOLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED.'
     .                    ,NCHAR,NERR,LEVEL,1,JP,IDUM,0,IDUM,IDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 190
          END IF
   50     CONTINUE
C     END PROCEDURE
          IF (CHECKL) THEN
C     DO(CHECK LENGTHS OF ARRAYS)
C     PROCEDURE(CHECK LENGTHS OF ARRAYS)
C
C     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE
C     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE.
              IF (LMDW.LT.MROWS) THEN
                  NERR = 11
                  NCHAR = 76
                  CALL XERRWV(
     .'DBOLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE.THE NUMBER OF
     . ROWS=(I2).',NCHAR,NERR,LEVEL,2,LMDW,MROWS,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LNDW.LT.NCOLS+1) THEN
                  NERR = 12
                  NCHAR = 69
                  CALL XERRWV(
     .'DBOLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+1=(
     .I2).',NCHAR,NERR,LEVEL,2,LNDW,NCOLS+1,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LLB.LT.NCOLS) THEN
                  NERR = 13
                  NCHAR = 88
                  CALL XERRWV(
     .'DBOLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) M
     .UST BE .GE. NCOLS=(I2).',NCHAR,NERR,LEVEL,2,LLB,NCOLS,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LLX.LT.LENX) THEN
                  NERR = 14
                  NCHAR = 70
                  CALL XERRWV(
     .'DBOLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. THE REQD. LENGTH=
     .(I2).',NCHAR,NERR,LEVEL,2,LLX,LENX,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LLRW.LT.5*NCOLS) THEN
                  NERR = 15
                  NCHAR = 62
                  CALL XERRWV(
     .  'DBOLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 5*NCOLS=(I2).'
     .                        ,NCHAR,NERR,LEVEL,2,LLRW,5*NCOLS,0,RDUM,
     .                        RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LLIW.LT.2*NCOLS) THEN
                  NERR = 16
                  NCHAR = 61
                  CALL XERRWV(
     .   'DBOLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS=(I2).'
     .                        ,NCHAR,NERR,LEVEL,2,LLIW,2*NCOLS,0,RDUM,
     .                        RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
              IF (LIOPT.LT.LP+1) THEN
                  NERR = 17
                  NCHAR = 71
                  CALL XERRWV(
     .'DBOLS(). THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. LEN.
     .=(I2).',NCHAR,NERR,LEVEL,2,LIOPT,LP+1,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 190
              END IF
C     END PROCEDURE
          END IF
      END IF
      GO TO (60,90),IGO
      GO TO 180
C
C     GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES
C     EQUATIONS AND DIRECTIONS TO QUIT PROCESSING.
C     CASE 1
   60 CONTINUE
C     DO(ACCUMULATE LEAST SQUARES EQUATIONS)
C     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS)
      MROWS = IOPT(LOCACC+1) - 1
      INROWS = IOPT(LOCACC+2)
      MNEW = MROWS + INROWS
      IF (MNEW.LT.0 .OR. MNEW.GT.MDW) THEN
          NERR = 10
          NCHAR = 61
          CALL XERRWV(
     .   'DBOLS(). NO. OF ROWS=(I1) MUST BE .GE. 0 .AND. .LE. MDW=(I2).'
     .                ,NCHAR,NERR,LEVEL,2,MNEW,MDW,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 190
      END IF
      DO 80 J = 1,MIN(NCOLS+1,MNEW)
          DO 70 I = MNEW,MAX(MROWS,J) + 1,-1
              IBIG = IDAMAX(I-J,W(J,J),1) + J - 1
C
C     PIVOT FOR INCREASED STABILITY.
              CALL DROTG(W(IBIG,J),W(I,J),SC,SS)
              CALL DROT(NCOLS+1-J,W(IBIG,J+1),MDW,W(I,J+1),MDW,SC,SS)
              W(I,J) = ZERO
   70     CONTINUE
   80 CONTINUE
      MROWS = MIN(NCOLS+1,MNEW)
      IOPT(LOCACC+1) = MROWS + 1
      IGO = IOPT(LOCACC)
C     END PROCEDURE
      IF (IGO.EQ.2) THEN
          IGO = 0
      END IF
      GO TO 180
C     CASE 2
   90 CONTINUE
C     DO(INITIALIZE VARIABLES AND DATA VALUES)
C     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES)
      DO 150 J = 1,NCOLS
          GO TO (100,110,120,130),ISCALE
          GO TO 140
  100     CONTINUE
C     CASE 1
C
C     THIS IS THE NOMINAL SCALING. EACH NONZERO
C     COL. HAS MAX. NORM EQUAL TO ONE.
          IBIG = IDAMAX(MROWS,W(1,J),1)
          RW(J) = ABS(W(IBIG,J))
          IF (RW(J).EQ.ZERO) THEN
              RW(J) = ONE
          ELSE
              RW(J) = ONE/RW(J)
          END IF
          GO TO 140
  110     CONTINUE
C     CASE 2
C
C     THIS CHOICE OF SCALING MAKES EACH NONZERO COLUMN
C     HAVE EUCLIDEAN LENGTH EQUAL TO ONE.
          RW(J) = DNRM2(MROWS,W(1,J),1)
          IF (RW(J).EQ.ZERO) THEN
              RW(J) = ONE
          ELSE
              RW(J) = ONE/RW(J)
          END IF
          GO TO 140
  120     CONTINUE
C     CASE 3
C
C     THIS CASE EFFECTIVELY SUPPRESSES SCALING BY SETTING
C     THE SCALING MATRIX TO THE IDENTITY MATRIX.
          RW(1) = ONE
          CALL DCOPY(NCOLS,RW,0,RW,1)
          GO TO 160
  130     CONTINUE
C     CASE 4
          GO TO 160
  140     CONTINUE
  150 CONTINUE
  160 CONTINUE
C     END PROCEDURE
C     DO(SOLVE BOUNDED LEAST SQUARES PROBLEM)
C     PROCEDURE(SOLVE BOUNDED LEAST SQUARES PROBLEM)
C
C     INITIALIZE IBASIS(*), J=1,NCOLS, AND IBB(*), J=1,NCOLS,
C     TO =J,AND =1, FOR USE IN DBOLSM( ).
      DO 170 J = 1,NCOLS
          IW(J) = J
          IW(J+NCOLS) = 1
          RW(3*NCOLS+J) = BL(J)
          RW(4*NCOLS+J) = BU(J)
  170 CONTINUE
      CALL DBOLSM(W,MDW,MROWS,NCOLS,RW(3*NCOLS+1),RW(4*NCOLS+1),IND,
     .            IOPT(LOPT),X,RNORM,MODE,RW(NCOLS+1),RW(2*NCOLS+1),RW,
     .            IW,IW(NCOLS+1))
C     END PROCEDURE
      IGO = 0
  180 CONTINUE
      RETURN
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  190 IF(MODE.GE.0)MODE = -NERR
      IGO = 0
      RETURN
C     END PROCEDURE
      END
      SUBROUTINE DBOLSM(W,MDW,MINPUT,NCOLS,BL,BU,IND,IOPT,X,RNORM,MODE,
     .                  RW,WW,SCL,IBASIS,IBB)
C***BEGIN PROLOGUE  DBOLSM
C***REFER TO  DBOCLS,DBOLS
C***ROUTINES CALLED  IVOUT,D1MACH,DAXPY,DCOPY,DDOT,DMOUT,DNRM2,DROT,
C                    DROTG,DSWAP,DVOUT,XERRWV
C***COMMON BLOCKS    DBOCOM
C***DESCRIPTION
C
C          Solve E*X = F (least squares sense) with bounds on
C            selected X values.
C     The user must have dimension statements of the form:
C
C       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS),
C      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS)
C       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS)
C
C     (here NX=number of extra locations required for options 1,...,7;
C     NX=0 for no options; here NI=number of extra locations possibly
C     required for options 1-7; NI=0 for no options; NI=14 if all the
C     options are simultaneously in use.)
C
C    INPUT
C    -----
C
C    --------------------
C    W(MDW,*),MROWS,NCOLS
C    --------------------
C     The array w(*,*) contains the matrix [E:F] on entry. The matrix
C     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in
C     the array W(*,*) with E occupying the first NCOLS columns and the
C     right side vector F in column NCOLS+1. The row dimension, MDW, of
C     the array W(*,*) must satisfy the inequality MDW .ge. MROWS.
C     Other values of MDW are errors. The values of MROWS and NCOLS
C     must be positive. Other values are errors.
C
C    ------------------
C    BL(*),BU(*),IND(*)
C    ------------------
C     These arrays contain the information about the bounds that the
C     solution values are to satisfy. The value of IND(J) tells the
C     type of bound and BL(J) and BU(J) give the explicit values for
C     the respective upper and lower bounds.
C
C    1.    For IND(J)=1, require X(J) .ge. BL(J).
C    2.    For IND(J)=2, require X(J) .le. BU(J).
C    3.    For IND(J)=3, require X(J) .ge. BL(J) and
C                                X(J) .le. BU(J).
C    4.    For IND(J)=4, no bounds on X(J) are required.
C     The values of BL(*),BL(*) are modified by the subprogram. Values
C     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3
C     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an
C     error.
C
C    -------
C    IOPT(*)
C    -------
C     This is the array where the user can specify nonstandard options
C     for DBOLSM( ). Most of the time this feature can be ignored by
C     setting the input value IOPT(1)=99. Occasionally users may have
C     needs that require use of the following subprogram options. For
C     details about how to use the options see below: IOPT(*) CONTENTS.
C
C     Option Number   Brief Statement of Purpose
C     ----- ------   ----- --------- -- -------
C           1         Move the IOPT(*) processing pointer.
C           2         Change rank determination tolerance.
C           3         Change blow-up factor that determines the
C                     size of variables being dropped from active
C                     status.
C           4         Reset the maximum number of iterations to use
C                     in solving the problem.
C           5         The data matrix is triangularized before the
C                     problem is solved whenever (NCOLS/MROWS) .lt.
C                     FAC. Change the value of FAC.
C           6         Redefine the weighting matrix used for
C                     linear independence checking.
C           7         Debug output is desired.
C          99         No more options to change.
C
C    ----
C    X(*)
C    ----
C     This array is used to pass data associated with options 1,2,3 and
C     5. Ignore this input parameter if none of these options are used.
C     Otherwise see below: IOPT(*) CONTENTS.
C
C    ----------------
C    IBASIS(*),IBB(*)
C    ----------------
C     These arrays must be initialized by the user. The values
C         IBASIS(J)=J, J=1,...,NCOLS
C         IBB(J)   =1, J=1,...,NCOLS
C     are appropriate except when using nonstandard features.
C
C    ------
C    SCL(*)
C    ------
C     This is the array of scaling factors to use on the columns of the
C     matrix E. These values must be defined by the user. To suppress
C     any column scaling set SCL(J)=1.0, J=1,...,NCOLS.
C
C    OUTPUT
C    ------
C
C    ----------
C    X(*),RNORM
C    ----------
C     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for
C     the constrained least squares problem. The value RNORM is the
C     minimum residual vector length.
C
C    ----
C    MODE
C    ----
C     The sign of mode determines whether the subprogram has completed
C     normally, or encountered an error condition or abnormal status.
C     A value of MODE .ge. 0 signifies that the subprogram has completed
C     normally. The value of MODE (.ge. 0) is the number of variables
C     in an active status: not at a bound nor at the value ZERO, for
C     the case of free variables. A negative value of MODE will be one
C     of the 20 cases -40,-39,...,-22, or -1. Values .lt. -1 correspond
C     to an abnormal completion of the subprogram. To understand the
C     abnormal completion codes see below: ERROR MESSAGES for DBOLSM( )
C     An approximate solution will be returned to the user only when
C     max. iterations is reached, MODE=-22.
C
C    -----------
C    RW(*),WW(*)
C    -----------
C     These are working arrays each with NCOLS entries. The array RW(*)
C     contains the working (scaled, nonactive) solution values. The
C     array WW(*) contains the working (scaled, active) gradient vector
C     values.
C
C    ----------------
C    IBASIS(*),IBB(*)
C    ----------------
C     These arrays contain information about the status of the solution
C     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the
C     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are
C     the active variables. The value (IBB(J)-1) is the number of times
C     variable J was reflected from its upper bound. (normally the user
C     can ignore these parameters.)
C
C    IOPT(*) CONTENTS
C    ------- --------
C     The option array allows a user to modify internal variables in
C     the subprogram without recompiling the source code. A central
C     goal of the initial software design was to do a good job for most
C     people. Thus the use of options will be restricted to a select
C     group of users. The processing of the option array proceeds as
C     follows: a pointer, here called LP, is initially set to the value
C     1. The value is updated as the options are processed.  At the
C     pointer position the option number is extracted and used for
C     locating other information that allows for options to be changed.
C     The portion of the array IOPT(*) that is used for each option is
C     fixed; the user and the subprogram both know how many locations
C     are needed for each option. A great deal of error checking is
C     done by the subprogram on the contents of the option array.
C     Nevertheless it is still possible to give the subprogram optional
C     input that is meaningless. For example some of the options use
C     the location X(NCOLS+IOFF) for passing data. The user must manage
C     the allocation of these locations when more than one piece of
C     option data is being passed to the subprogram.
C
C   1
C   -
C     Move the processing pointer (either forward or backward) to the
C     location IOPT(LP+1). The processing pointer is moved to location
C     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over
C     locations 3,...,NCOLS+2 of IOPT(*),
C
C       IOPT(1)=1
C       IOPT(2)=NCOLS+3
C       (IOPT(I), I=3,...,NCOLS+2 are not defined here.)
C       IOPT(NCOLS+3)=99
C       CALL DBOLSM( )
C
C     CAUTION: Misuse of this option can yield some very hard
C     -to-find bugs.  Use it with care.
C
C   2
C   -
C     The algorithm that solves the bounded least squares problem
C     iteratively drops columns from the active set. This has the
C     effect of joining a new column vector to the QR factorization of
C     the rectangular matrix consisting of the partially triangularized
C     nonactive columns. After triangularizing this matrix a test is
C     made on the size of the pivot element. The column vector is
C     rejected as dependent if the magnitude of the pivot element is
C     .le. TOL* magnitude of the column in components strictly above
C     the pivot element. Nominally the value of this (rank) tolerance
C     is TOL = DRELPR, where DRELPR is relative machine
C     precision. To change only the value of TOL, for example,
C
C       X(NCOLS+1)=TOL
C       IOPT(1)=2
C       IOPT(2)=1
C       IOPT(3)=99
C       CALL DBOLSM()
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=TOL
C       IOPT(LP)=2
C       IOPT(LP+1)=IOFF
C        .
C       CALL DBOLSM()
C
C     The required length of IOPT(*) is increased by 2 if option 2 is
C     used; The required length of X(*) is increased by 1. A value of
C     IOFF .le. 0 is an error. A value of TOL .le. DRELPR gives a
C     warning message; it is not considered an error.
C     Here DRELPR is the relative machine precision.
C
C   3
C   -
C     A solution component is left active (not used) if, roughly
C     speaking, it seems too large. Mathematically the new component is
C     left active if the magnitude is .ge.((vector norm of F)/(matrix
C     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(DRELPR)
C     where DRELPR is the relative machine precision. To change only
C     the value of BLOWUP, for example,
C
C       X(NCOLS+2)=BLOWUP
C       IOPT(1)=3
C       IOPT(2)=2
C       IOPT(3)=99
C       CALL DBOLSM()
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=BLOWUP
C       IOPT(LP)=3
C       IOPT(LP+1)=IOFF
C        .
C       CALL DBOLSM()
C
C     The required length of IOPT(*) is increased by 2 if option 3 is
C     used; the required length of X(*) is increased by 1. A value of
C     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error.
C
C   4
C   -
C     Normally the algorithm for solving the bounded least squares
C     problem requires between NCOLS/3 and NCOLS drop-add steps to
C     converge. (this remark is based on examining a small number of
C     test cases.) The amount of arithmetic for such problems is
C     typically about twice that required for linear least squares if
C     there are no bounds and if plane rotations are used in the
C     solution method. Convergence of the algorithm, while
C     mathematically certain, can be much slower than indicated. To
C     avoid this potential but unlikely event ITMAX drop-add steps are
C     permitted. Nominally ITMAX=5*(MAX(MROWS,NCOLS)). To change the
C     value of ITMAX, for example,
C
C       IOPT(1)=4
C       IOPT(2)=ITMAX
C       IOPT(3)=99
C       CALL DBOLSM()
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       IOPT(LP)=4
C       IOPT(LP+1)=ITMAX
C        .
C       CALL DBOLSM()
C
C     The value of ITMAX must be .gt. 0. Other values are errors. Use
C     of this option increases the required length of IOPT(*) by 2.
C
C   5
C   -
C     For purposes of increased efficiency the MROWS by NCOLS+1 data
C     matrix [E:F] is triangularized as a first step whenever MROWS
C     satisfies FAC*MROWS .gt. NCOLS. Nominally FAC=0.  To change the
C     value of FAC,
C
C       X(NCOLS+3)=FAC
C       IOPT(1)=5
C       IOPT(2)=3
C       IOPT(3)=99
C       CALL DBOLSM()
C
C     Generally, if LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=FAC
C       IOPT(LP)=5
C       IOPT(LP+1)=IOFF
C        .
C       CALL DBOLSM()
C
C     The value of FAC must be nonnegative. Other values are errors.
C     Using the value FAC=0.0 suppresses the initial triangularization step.
C     Use of this option increases the required length of IOPT(*) by 2;
C     The required length of of X(*) is increased by 1.
C
C   6
C   -
C     The norm used in testing the magnitudes of the pivot element
C     compared to the mass of the column above the pivot line can be
C     changed. The type of change that this option allows is to weight
C     the components with an index larger than MVAL by the parameter
C     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and
C     WT, where LP is the processing pointer for IOPT(*),
C
C       X(NCOLS+IOFF)=WT
C       IOPT(LP)=6
C       IOPT(LP+1)=IOFF
C       IOPT(LP+2)=MVAL
C
C     Use of this option increases the required length of IOPT(*) by 3.
C     The length of X(*) is increased by 1. Values of MVAL must be
C     nonnegative and not greater than MROWS. Other values are errors.
C     The value of WT must be positive. Any other value is an error. If
C     either error condition is present a message will be printed.
C
C   7
C   -
C     Debug output, showing the detailed add-drop steps for the
C     constrained least squares problem, is desired. This option is
C     intended to be used to locate suspected bugs.  To print,
C
C       IOPT(LP)=7
C
C   99
C   --
C     There are no more options to change.
C
C     The values for options are 1,...,7,99, and are the only ones
C     permitted. Other values are errors. Options -99,-1,...,-7 mean
C     that the repective options 99,1,...,7 are left at their default
C     values. An example is the option to modify the (rank) tolerance:
C
C       X(NCOLS+1)=TOL
C       IOPT(1)=-2
C       IOPT(2)=1
C       IOPT(3)=99
C
C    Error Messages for DBOLSM( )
C    ----- -------- --- ---------
C    WARNING IN...
C    DBOLSM(). MORE THAN (I1)=ITMAX ITERATIONS SOLVING BOUNDED LEAST
C    SQUARES PROBLEM.
C              IN ABOVE MESSAGE, I1=         3
C    ERROR NUMBER =        22
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM. THE OPTION NUMBER=(I1) IS NOT DEFINED.
C              IN ABOVE MESSAGE, I1=         0
C    ERROR NUMBER =        23
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE OFFSET=(I1) BEYOND POSTION NCOLS=(I2)
C    MUST BE POSITIVE FOR OPTION NUMBER 2.
C              IN ABOVE MESSAGE, I1=         0
C              IN ABOVE MESSAGE, I2=         1
C    ERROR NUMBER =        24
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE TOLERANCE FOR RANK DETERMINATION=(R1)
C    IS LESS THAN MACHINE PRECISION=(R2).
C              IN ABOVE MESSAGE, R1=    0.
C              IN ABOVE MESSAGE, R2=     .7105427358E-14
C    ERROR NUMBER =        25
C
C    WARNING IN...
C    DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST
C     BE POSTIVE FOR OPTION NUMBER 3.
C              IN ABOVE MESSAGE, I1=         0
C              IN ABOVE MESSAGE, I2=         1
C    ERROR NUMBER =        26
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING
C    VARIABLES MUST BE POSITIVE. NOW=(R1).
C              IN ABOVE MESSAGE, R1=    0.
C    ERROR NUMBER =        27
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE MAXIMUM NUMBER OF ITERATIONS=(I1) MUST BE POSITIVE.
C              IN ABOVE MESSAGE, I1=         0
C    ERROR NUMBER =        28
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST BE
C    POSTIVE FOR OPTION NUMBER 5.
C              IN ABOVE MESSAGE, I1=         0
C              IN ABOVE MESSAGE, I2=         1
C    ERROR NUMBER =        29
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE FACTOR (NCOLS/MROWS) WHERE PRETRIANGULARIZING IS
C    PERFORMED MUST BE NONNEGATIVE. NOW=(R1).
C              IN ABOVE MESSAGE, R1=    -.2500000000E+00
C    ERROR NUMBER =        30
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE NUMBER OF ROWS=(I1) MUST BE POSITIVE.
C              IN ABOVE MESSAGE, I1=         0
C    ERROR NUMBER =        31
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE NUMBER OF COLS.=(I1) MUST BE POSTIVE.
C              IN ABOVE MESSAGE, I1=         0
C    ERROR NUMBER =        32
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE
C    NUMBER OF ROWS =(I2).
C              IN ABOVE MESSAGE, I1=         0
C              IN ABOVE MESSAGE, I2=         1
C    ERROR NUMBER =        33
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). FOR J=(I1) THE CONSTRAINT INDICATOR MUST BE 1-4.
C              IN ABOVE MESSAGE, I1=         1
C              IN ABOVE MESSAGE, I2=         0
C    ERROR NUMBER =        34
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). FOR J=(I1) THE LOWER BOUND=(R1) IS .GT. THE UPPER
C     BOUND=(R2).
C              IN ABOVE MESSAGE, I1=         1
C              IN ABOVE MESSAGE, R1=    0.
C              IN ABOVE MESSAGE, R2=    -.1000000000E+01
C    ERROR NUMBER =        35
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE INPUT ORDER OF COLUMNS=(I1) IS NOT BETWEEN 1
C    AND NCOLS=(I2).
C              IN ABOVE MESSAGE, I1=         0
C              IN ABOVE MESSAGE, I2=         1
C    ERROR NUMBER =        36
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE BOUND POLARITY FLAG IN COMPONENT J=(I1) MUST
C    BE POSITIVE.  NOW=(I2).
C              IN ABOVE MESSAGE, I1=         1
C              IN ABOVE MESSAGE, I2=         0
C    ERROR NUMBER =        37
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE ROW SEPARATOR TO APPLY WEIGHTING (I1) MUST LIE
C    BETWEEN 0 AND MROWS (I2). WEIGHT (R1) MUST BE POSITIVE.
C              IN ABOVE MESSAGE, I1=        -1
C              IN ABOVE MESSAGE, I2=         2
C              IN ABOVE MESSAGE, R1=    0.
C    ERROR NUMBER =        38
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE OFFSET (I1) BEYOND POSITION NCOLS=(I2) MUST BE
C    POSITIVE FOR OPTION NUMBER 7.
C              IN ABOVE MESSAGE, I1=        -1
C              IN ABOVE MESSAGE, I2=         2
C    ERROR NUMBER =        39
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C
C    WARNING IN...
C    DBOLSM(). THE COLUMN PIVOTING THRESHOLD FACTOR MUST BE
C    POSITIVE. NOW=(R1).
C              IN ABOVE MESSAGE, R1=    0.
C    ERROR NUMBER =        40
C    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.)
C***END PROLOGUE  DBOLSM
C
C     PURPOSE
C     -------
C     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED
C     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS:
C
C     SOLVE E*X =  F  (LEAST SQUARES SENSE)
C     WITH BOUNDS ON SELECTED X VALUES.
C
C     REVISED 880722-1100
C     REVISED YYMMDD-HHMM
C
C     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN
C     EDITING AT THE CARD 'C++'.
C     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS
C     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/,
C     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/,
C     /SROT/ TO /DROT/, /SROTG/ TO /DROTG/, /R1MACH/ TO /D1MACH/,
C     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/,
C     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/,
C     /REAL            / TO /DOUBLE PRECISION/.
C ++
C
C     THIS VARIABLE REMAINS TYPE REAL.
      REAL RDUM
      DOUBLE PRECISION W(MDW,*),BL(*),BU(*)
      DOUBLE PRECISION X(*),RW(*),WW(*),SCL(*)
      DOUBLE PRECISION ALPHA,BETA,BOU,COLABV,COLBLO
      DOUBLE PRECISION CL1,CL2,CL3,ONE,BIG
      DOUBLE PRECISION FAC,RNORM,SC,SS,T,TOLIND,WT
      DOUBLE PRECISION TWO,T1,T2,WLARGE,WLA,WLB,XNEW
      DOUBLE PRECISION ZERO,DDOT,DNRM2
      DOUBLE PRECISION D1MACH,TOLSZE
      INTEGER IBASIS(*),IBB(*),IND(*),IOPT(*)
      LOGICAL FOUND,CONSTR,CNZ
      COMMON /DBOCOM/IDOPE(5)
      INEXT(IDUM) = MIN(IDUM+1,MROWS)
C***FIRST EXECUTABLE STATEMENT  DBOLSM
      LEVEL = 1
C
C    VERIFY THAT THE PROBLEM DIMENSIONS ARE DEFINED PROPERLY.
      IF (MINPUT.LE.0) THEN
          NERR = 31
          NCHAR = 51
          CALL XERRWV(
     .             'DBOLSM(). THE NUMBER OF ROWS=(I1) MUST BE POSITIVE.'
     .                ,NCHAR,NERR,LEVEL,1,MINPUT,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 600

      END IF

      IF (NCOLS.LE.0) THEN
          NERR = 32
          NCHAR = 51
          CALL XERRWV(
     .             'DBOLSM(). THE NUMBER OF COLS.=(I1) MUST BE POSTIVE.'
     .                ,NCHAR,NERR,LEVEL,1,NCOLS,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 600

      END IF

      IF (MDW.LT.MINPUT) THEN
          NERR = 33
          NCHAR = 78
          CALL XERRWV(
     .'DBOLSM(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE NUMBER
     .OF ROWS=(I2).',NCHAR,NERR,LEVEL,2,MDW,MROWS,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 600

      END IF
C
C     VERIFY THAT BOUND INFORMATION IS CORRECT.
      DO 10 J = 1,NCOLS
         IF (IND(J).LT.1 .OR. IND(J).GT.4) THEN
             NERR = 34
             NCHAR = 58
             CALL XERRWV(
     .      'DBOLSM(). FOR J=(I1) THE CONSTRAINT INDICATOR MUST BE 1-4.'
     .                   ,NCHAR,NERR,LEVEL,2,J,IND(J),0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
             GO TO 600

         END IF

   10 CONTINUE
      DO 20 J = 1,NCOLS
         IF (IND(J).EQ.3) THEN
             IF (BU(J).LT.BL(J)) THEN
                 NERR = 35
                 NCHAR = 71
                 CALL XERRWV(
     .'DBOLSM(). FOR J=(I1) THE LOWER BOUND=(R1) IS .GT. THE UPPER BOUND
     .=(R2).',NCHAR,NERR,LEVEL,1,J,IDUM,2,BL(J),BU(J))
C     DO(RETURN TO USER PROGRAM UNIT)
                 GO TO 600

             END IF

         END IF
*
   20 CONTINUE
C
C     CHECK THAT PERMUTATION AND POLARITY ARRAYS HAVE BEEN SET.
      DO 30 J = 1,NCOLS
         IF (IBASIS(J).LT.1 .OR. IBASIS(J).GT.NCOLS) THEN
             NERR = 36
             NCHAR = 74
             CALL XERRWV(
     .'DBOLSM(). THE INPUT ORDER OF COLUMNS=(I1) IS NOT BETWEEN 1 AND NC
     .OLS=(I2).',NCHAR,NERR,LEVEL,2,IBASIS(J),NCOLS,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
             GO TO 600

         END IF

         IF (IBB(J).LE.0) THEN
             NERR = 37
             NCHAR = 81
             CALL XERRWV(
     .'DBOLSM(). THE BOUND POLARITY FLAG IN COMPONENT J=(I1) MUST BE POS
     .ITIVE. NOW=(I2).',NCHAR,NERR,LEVEL,2,J,IBB(J),0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
             GO TO 600

         END IF

   30 CONTINUE
C     DO(PROCESS OPTION ARRAY)
      GO TO 570

   40 CONTINUE
C     DO(INITIALIZE VARIABLES AND DATA VALUES)
      GO TO 460

   50 CONTINUE
      IF (IPRINT.GT.0) THEN
          CALL DMOUT(MROWS,NCOLS+1,MDW,W,'('' PRETRI. INPUT MATRIX'')',
     .               -4)
          CALL DVOUT(NCOLS,BL,'('' LOWER BOUNDS'')',-4)
          CALL DVOUT(NCOLS,BU,'('' UPPER BOUNDS'')',-4)
      END IF

   60 CONTINUE
      ITER = ITER + 1
      IF (ITER.LE.ITMAX) GO TO 80
      NERR = 22
      NCHAR = 80
      CALL XERRWV(
     .'DBOLSM(). MORE THAN (I1)=ITMAX ITERATIONS SOLVING BOUNDED LEAST S
     .QUARES PROBLEM.',NCHAR,NERR,LEVEL,1,ITMAX,IDUM,0,RDUM,RDUM)
C     DO(RESCALE AND TRANSLATE VARIABLES)
      IGOPR = 1
      GO TO 130

   70 CONTINUE
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 600

   80 CONTINUE
C     DO(FIND A VARIABLE TO BECOME NON-ACTIVE)
      GO TO 180

   90 CONTINUE
      IF (FOUND) GO TO 110
C     DO(RESCALE AND TRANSLATE VARIABLES)
      IGOPR = 2
      GO TO 130

  100 CONTINUE
      MODE = NSETB
      RETURN

  110 CONTINUE
C     DO(MAKE MOVE AND UPDATE FACTORIZATION)
      GO TO 280

  120 CONTINUE
      GO TO 60
C     PROCEDURE(RESCALE AND TRANSLATE VARIABLES)
  130 CONTINUE
      CALL DCOPY(NSETB,X,1,RW,1)
      X(1) = ZERO
      CALL DCOPY(NCOLS,X,0,X,1)
      DO 140 J = 1,NSETB
         JCOL = ABS(IBASIS(J))
         X(JCOL) = RW(J)*ABS(SCL(JCOL))
  140 CONTINUE
      DO 150 J = 1,NCOLS
         IF (MOD(IBB(J),2).EQ.0) X(J) = BU(J) - X(J)
  150 CONTINUE
      DO 160 J = 1,NCOLS
         JCOL = IBASIS(J)
         IF (JCOL.LT.0) X(-JCOL) = BL(-JCOL) + X(-JCOL)
  160 CONTINUE
      DO 170 J = 1,NCOLS
         IF (SCL(J).LT.ZERO) X(J) = -X(J)
  170 CONTINUE
      CALL DSCAL(MROWS-MVAL,WT,W(INEXT(MVAL),NCOLS+1),1)
      RNORM = DNRM2(MROWS-MAX(NSETB,MVAL),
     .        W(INEXT(MAX(NSETB,MVAL)),NCOLS+1),1)
C     END PROCEDURE
      GO TO (70,100),IGOPR
C     PROCEDURE(FIND A VARIABLE TO BECOME NON-ACTIVE)
  180 CONTINUE
C
C     COMPUTE (NEGATIVE) OF GRADIENT VECTOR, W=
C     (TRANSPOSE OF E)*(F-E*X).
      WW(1) = ZERO
      CALL DCOPY(NCOLS,WW,0,WW,1)
      DO 190 J = NSETB + 1,NCOLS
         JCOL = ABS(IBASIS(J))
         WW(J) = DDOT(MROWS-NSETB,W(INEXT(NSETB),J),1,
     .           W(INEXT(NSETB),NCOLS+1),1)*ABS(SCL(JCOL))
  190 CONTINUE
      IF (IPRINT.GT.0) THEN
          CALL DVOUT(NCOLS,WW,'('' GRADIENT VALUES'')',-4)
          CALL IVOUT(NCOLS,IBASIS,'('' INTERNAL VARIABLE ORDER'')',-4)
          CALL IVOUT(NCOLS,IBB,'('' BOUND POLARITY'')',-4)
      END IF

  200 CONTINUE
C
C     IF ACTIVE SET = NUMBER OF TOTAL ROWS, QUIT.
      IF (NSETB.EQ.MROWS) THEN
          FOUND = .FALSE.
C     EXIT PROCEDURE
          GO TO 90

      END IF
C
C     CHOOSE AN EXTREMAL COMPONENT OF GRADIENT VECTOR
C     FOR A CANDIDATE TO BECOME NON-ACTIVE.
      WLARGE = -BIG
      JBIG = 0
      CNZ = .FALSE.
      DO 210 J = NSETB + 1,NCOLS
         T = WW(J)
C     SKIP LOOKING AT COMPONENTS FLAGGED AS NON-CANDIDATES.
         IF (T.EQ.BIG) GO TO 210
         ITEMP = IBASIS(J)
         JCOL = ABS(ITEMP)
         IF (NSETB.LT.MVAL) THEN
             CL1 = DNRM2(NSETB,W(1,J),1)
             CL2 = DNRM2(MVAL-NSETB,W(INEXT(NSETB),J),1)
             COLABV = CL1
             COLBLO = CL2

         ELSE
             CL1 = DNRM2(MVAL,W(1,J),1)
             CL2 = ABS(WT)*DNRM2(NSETB-MVAL,W(INEXT(MVAL),J),1)
             CL3 = ABS(WT)*DNRM2(MROWS-NSETB,W(INEXT(NSETB),J),1)
             CALL DROTG(CL1,CL2,SC,SS)
             COLABV = ABS(CL1)
             COLBLO = CL3
         END IF

         IF (ITEMP.LT.0) THEN
             IF (MOD(IBB(JCOL),2).EQ.0) T = -T
C     SKIP LOOKING AT COMPONENTS THAT WOULD NOT DECREASE OBJECTIVE.
             IF (T.LT.ZERO) GO TO 210
         END IF
C     THIS IS A COLUMN PIVOTING STEP THAT MAXIMIZES THE RATIO OF
C     COLUMN MASS ON AND BELOW THE PIVOT LINE RELATIVE TO THAT
C     STRICTLY ABOVE THE PIVOT LINE.
         IF (COLABV.EQ.ZERO .AND. .NOT. CNZ) THEN
             T = COLBLO*ABS(SCL(JCOL))
             IF (WLARGE.LT.T) THEN
                 WLARGE = T
                 JBIG = J
             END IF

         ELSE
             IF ( .NOT. CNZ) THEN
                 WLA = ZERO
                 WLB = ZERO
                 CNZ = .TRUE.
             END IF

           IF(SQRT(COLBLO)*SQRT(WLA) .GE. SQRT(COLABV)*SQRT(WLB)) THEN
                WLB=COLBLO
                WLA=COLABV
                JBIG=J
           END IF
         END IF

  210 CONTINUE
      IF (JBIG.EQ.0) THEN
          FOUND = .FALSE.
          IF (IPRINT.GT.0) THEN
              CALL IVOUT(0,I,'('' FOUND NO VARIABLE TO ENTER'')',-4)
          END IF
C     EXIT PROCEDURE
          GO TO 90

      END IF
C
C     SEE IF THE INCOMING COL. IS SUFFICIENTLY INDEPENDENT.
C     THIS TEST IS MADE BEFORE AN ELIMINATION IS PERFORMED.
      IF (IPRINT.GT.0) THEN
          CALL IVOUT(1,JBIG,'('' TRY TO BRING IN THIS COL.'')',-4)
      END IF

      IF (CNZ) THEN
      IF(WLB.LE.WLA*TOLIND) THEN
          FOUND = .FALSE.
          IF (IPRINT.GT.0) THEN
              CALL IVOUT(0,I,'('' VARIABLE IS DEPENDENT, NOT USED.'')',
     .                   -4)
          END IF

          WW(JBIG) = BIG
          GO TO 200

      END IF
      END IF
C
C     SWAP MATRIX COLS. NSETB+1 AND JBIG, PLUS POINTER INFO., AND
C     GRADIENT VALUES.
      NSETB = NSETB + 1
      IF (NSETB.NE.JBIG) THEN
          CALL DSWAP(MROWS,W(1,NSETB),1,W(1,JBIG),1)
          CALL DSWAP(1,WW(NSETB),1,WW(JBIG),1)
          ITEMP = IBASIS(NSETB)
          IBASIS(NSETB) = IBASIS(JBIG)
          IBASIS(JBIG) = ITEMP
      END IF
C
C     ELIMINATE ENTRIES BELOW THE PIVOT LINE IN COL. NSETB.
      IF (MROWS.GT.NSETB) THEN
          DO 220 I = MROWS,NSETB + 1,-1
             IF (I.EQ.MVAL+1) GO TO 220
             CALL DROTG(W(I-1,NSETB),W(I,NSETB),SC,SS)
             W(I,NSETB) = ZERO
             CALL DROT(NCOLS-NSETB+1,W(I-1,NSETB+1),MDW,W(I,NSETB+1),
     .                 MDW,SC,SS)
  220     CONTINUE
          IF ((MVAL.GE.NSETB) .AND. (MVAL.LT.MROWS)) THEN
              T = W(NSETB,NSETB)
              IF (T.NE.ZERO) THEN
                  T = WT*W(MVAL+1,NSETB)/T

              ELSE
                  T = BIG
              END IF

              IF (TOLIND*ABS(T).GT.ONE) THEN
                  CALL DSWAP(NCOLS-NSETB+2,W(NSETB,NSETB),MDW,
     .                       W(MVAL+1,NSETB),MDW)
                  CALL DSCAL(NCOLS-NSETB+2,WT,W(NSETB,NSETB),MDW)
                  CALL DSCAL(NCOLS-NSETB+2,ONE/WT,W(MVAL+1,NSETB),MDW)
              END IF

              CALL DROTG(W(NSETB,NSETB),W(MVAL+1,NSETB),SC,SS)
              W(MVAL+1,NSETB) = ZERO
              CALL DROT(NCOLS-NSETB+1,W(NSETB,NSETB+1),MDW,
     .                  W(MVAL+1,NSETB+1),MDW,SC,SS)
          END IF

      END IF

      IF (W(NSETB,NSETB).EQ.ZERO) THEN
          WW(NSETB) = BIG
          NSETB = NSETB - 1
          IF (IPRINT.GT.0) THEN
              CALL IVOUT(0,I,'('' PIVOT IS ZERO, NOT USED.'')',-4)
          END IF

          GO TO 200

      END IF
C
C     CHECK THAT NEW VARIABLE IS MOVING IN THE RIGHT DIRECTION.
      ITEMP = IBASIS(NSETB)
      JCOL = ABS(ITEMP)
      XNEW = (W(NSETB,NCOLS+1)/W(NSETB,NSETB))/ABS(SCL(JCOL))
C CONT: DO BLOCK
C QUIT: DO BLOCK
      IF (ITEMP.LT.0) THEN
C     IF(WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) EXIT(QUIT)
C     IF(WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO) EXIT(QUIT)
          IF (WW(NSETB).GE.ZERO .AND. XNEW.LE.ZERO) GO TO 230
          IF (WW(NSETB).LE.ZERO .AND. XNEW.GE.ZERO) GO TO 230
      END IF
C     EXIT(CONT)
      GO TO 240
C     END BLOCK
  230 CONTINUE
      WW(NSETB) = BIG
      NSETB = NSETB - 1
      IF (IPRINT.GT.0) THEN
          CALL IVOUT(0,I,'('' VARIABLE HAS BAD DIRECTION, NOT USED.'')',
     .               -4)
      END IF

      GO TO 200
C     END BLOCK
  240 CONTINUE
      FOUND = .TRUE.
C     EXIT PROCEDURE
      GO TO 250

  250 CONTINUE
C     END PROCEDURE
      GO TO 90
C     PROCEDURE(SOLVE THE TRIANGULAR SYSTEM)
  260 CONTINUE
      CALL DCOPY(NSETB,W(1,NCOLS+1),1,RW,1)
      DO 270 J = NSETB,1,-1
         RW(J) = RW(J)/W(J,J)
         JCOL = ABS(IBASIS(J))
         T = RW(J)
         IF (MOD(IBB(JCOL),2).EQ.0) RW(J) = -RW(J)
         CALL DAXPY(J-1,-T,W(1,J),1,RW,1)
         RW(J) = RW(J)/ABS(SCL(JCOL))
  270 CONTINUE
      IF (IPRINT.GT.0) THEN
          CALL DVOUT(NSETB,RW,'('' SOLN. VALUES'')',-4)
          CALL IVOUT(NSETB,IBASIS,'('' COLS. USED'')',-4)
      END IF
C     END PROCEDURE
      GO TO (290,430),LGOPR
C     PROCEDURE(MAKE MOVE AND UPDATE FACTORIZATION)
  280 CONTINUE
C     DO(SOLVE THE TRIANGULAR SYSTEM)
      LGOPR = 1
      GO TO 260

  290 CONTINUE
C
C     SEE IF THE UNCONSTRAINED SOL. (OBTAINED BY SOLVING THE
C     TRIANGULAR SYSTEM) SATISFIES THE PROBLEM BOUNDS.
      ALPHA = TWO
      BETA = TWO
      X(NSETB) = ZERO
      DO 300 J = 1,NSETB
         ITEMP = IBASIS(J)
         JCOL = ABS(ITEMP)
         T1 = TWO
         T2 = TWO
         IF (ITEMP.LT.0) THEN
             BOU = ZERO

         ELSE
             BOU = BL(JCOL)
         END IF

         IF ((-BOU).NE.BIG) BOU = BOU/ABS(SCL(JCOL))
         IF (RW(J).LE.BOU) T1 = (X(J)-BOU)/ (X(J)-RW(J))
         BOU = BU(JCOL)
         IF (BOU.NE.BIG) BOU = BOU/ABS(SCL(JCOL))
         IF (RW(J).GE.BOU) T2 = (BOU-X(J))/ (RW(J)-X(J))
C
C     IF NOT, THEN COMPUTE A STEP LENGTH SO THAT THE
C     VARIABLES REMAIN FEASIBLE.
         IF (T1.LT.ALPHA) THEN
             ALPHA = T1
             JDROP1 = J
         END IF

         IF (T2.LT.BETA) THEN
             BETA = T2
             JDROP2 = J
         END IF

  300 CONTINUE
      CONSTR = ALPHA .LT. TWO .OR. BETA .LT. TWO
      IF (CONSTR) GO TO 310
C
C     ACCEPT THE CANDIDATE BECAUSE IT SATISFIES THE STATED BOUNDS
C     ON THE VARIABLES.
      CALL DCOPY(NSETB,RW,1,X,1)
      GO TO 120

  310 CONTINUE
C
C     TAKE A STEP THAT IS AS LARGE AS POSSIBLE WITH ALL
C     VARIABLES REMAINING FEASIBLE.
      DO 320 J = 1,NSETB
         X(J) = X(J) + MIN(ALPHA,BETA)* (RW(J)-X(J))
  320 CONTINUE
      IF (ALPHA.LE.BETA) THEN
          JDROP2 = 0

      ELSE
          JDROP1 = 0
      END IF

  330 IF (JDROP1+JDROP2.GT.0 .AND. NSETB.GT.0) GO TO 340
      GO TO 450

  340 JDROP = JDROP1 + JDROP2
      ITEMP = IBASIS(JDROP)
      JCOL = ABS(ITEMP)
      IF (JDROP2.GT.0) THEN
C
C     VARIABLE IS AT AN UPPER BOUND.  SUBTRACT MULTIPLE OF THIS COL.
C     FROM RIGHT HAND SIDE.
          T = BU(JCOL)
          IF (ITEMP.GT.0) THEN
              BU(JCOL) = T - BL(JCOL)
              BL(JCOL) = -T
              ITEMP = -ITEMP
              SCL(JCOL) = -SCL(JCOL)
              DO 350 I = 1,JDROP
                 W(I,JDROP) = -W(I,JDROP)
  350         CONTINUE

          ELSE
              IBB(JCOL) = IBB(JCOL) + 1
              IF (MOD(IBB(JCOL),2).EQ.0) T = -T
          END IF
C     VARIABLE IS AT A LOWER BOUND.
      ELSE
          IF (ITEMP.LT.ZERO) THEN
              T = ZERO

          ELSE
              T = -BL(JCOL)
              BU(JCOL) = BU(JCOL) + T
              ITEMP = -ITEMP
          END IF

      END IF

      CALL DAXPY(JDROP,T,W(1,JDROP),1,W(1,NCOLS+1),1)
C
C     MOVE CERTAIN COLS. LEFT TO ACHIEVE UPPER HESSENBERG FORM.
      CALL DCOPY(JDROP,W(1,JDROP),1,RW,1)
      DO 360 J = JDROP + 1,NSETB
         IBASIS(J-1) = IBASIS(J)
         X(J-1) = X(J)
         CALL DCOPY(J,W(1,J),1,W(1,J-1),1)
  360 CONTINUE
      IBASIS(NSETB) = ITEMP
      W(1,NSETB) = ZERO
      CALL DCOPY(MROWS-JDROP,W(1,NSETB),0,W(JDROP+1,NSETB),1)
      CALL DCOPY(JDROP,RW,1,W(1,NSETB),1)
C
C     TRANSFORM THE MATRIX FROM UPPER HESSENBERG FORM TO
C     UPPER TRIANGULAR FORM.
      NSETB = NSETB - 1
C SMLL:
C     *DO BLOCK
C NRML:
C     *DO BLOCK
      DO 380 I = JDROP,NSETB
C
C     LOOK FOR SMALL PIVOTS AND AVOID MIXING WEIGHTED AND
C     NONWEIGHTED ROWS.
         IF (I.EQ.MVAL) THEN
             T = ZERO
             DO 370 J = I,NSETB
                JCOL = ABS(IBASIS(J))
                T1 = ABS(W(I,J)*SCL(JCOL))
                IF (T1.GT.T) THEN
                    JBIG = J
                    T = T1
                END IF

  370        CONTINUE
C     EXIT(NRML)
             GO TO 390

         END IF

         CALL DROTG(W(I,I),W(I+1,I),SC,SS)
         W(I+1,I) = ZERO
         CALL DROT(NCOLS-I+1,W(I,I+1),MDW,W(I+1,I+1),MDW,SC,SS)
  380 CONTINUE
C     EXIT (SMLL)
      GO TO 420
C     END BLOCK
  390 CONTINUE
C
C     THE TRIANGULARIZATION IS COMPLETED BY GIVING UP
C     THE HESSENBERG FORM AND TRIANGULARIZING A RECTANGULAR MATRIX.
      CALL DSWAP(MROWS,W(1,I),1,W(1,JBIG),1)
      CALL DSWAP(1,WW(I),1,WW(JBIG),1)
      CALL DSWAP(1,X(I),1,X(JBIG),1)
      ITEMP = IBASIS(I)
      IBASIS(I) = IBASIS(JBIG)
      IBASIS(JBIG) = ITEMP
      JBIG = I
      DO 410 J = JBIG,NSETB
         DO 400 I = J + 1,MROWS
            CALL DROTG(W(J,J),W(I,J),SC,SS)
            W(I,J) = ZERO
            CALL DROT(NCOLS-J+1,W(J,J+1),MDW,W(I,J+1),MDW,SC,SS)
  400    CONTINUE
  410 CONTINUE
C     END BLOCK
  420 CONTINUE
C
C     SEE IF THE REMAINING COEFFICIENTS ARE FEASIBLE.  THEY SHOULD
C     BE BECAUSE OF THE WAY MIN(ALPHA,BETA) WAS CHOSEN.  ANY THAT ARE
C     NOT FEASIBLE WILL BE SET TO THEIR BOUNDS AND
C     APPROPRIATELY TRANSLATED.
      JDROP1 = 0
      JDROP2 = 0
C     DO(SOLVE THE TRIANGULAR SYSTEM)
      LGOPR = 2
      GO TO 260

  430 CONTINUE
      CALL DCOPY(NSETB,RW,1,X,1)
      DO 440 J = 1,NSETB
         ITEMP = IBASIS(J)
         JCOL = ABS(ITEMP)
         IF (ITEMP.LT.0) THEN
             BOU = ZERO

         ELSE
             BOU = BL(JCOL)
         END IF

         IF ((-BOU).NE.BIG) BOU = BOU/ABS(SCL(JCOL))
         IF (X(J).LE.BOU) THEN
             JDROP1 = J
             GO TO 330

         END IF

         BOU = BU(JCOL)
         IF (BOU.NE.BIG) BOU = BOU/ABS(SCL(JCOL))
         IF (X(J).GE.BOU) THEN
             JDROP2 = J
             GO TO 330

         END IF

  440 CONTINUE
      GO TO 330

  450 CONTINUE
C     END PROCEDURE
      GO TO 120
C     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES)
  460 CONTINUE
C
C     PRETRIANGULARIZE RECTANGULAR ARRAYS OF CERTAIN SIZES
C     FOR INCREASED EFFICIENCY.
      IF (FAC*MINPUT.GT.NCOLS) THEN
          DO 480 J = 1,NCOLS + 1
             DO 470 I = MINPUT,J + MVAL + 1,-1
                CALL DROTG(W(I-1,J),W(I,J),SC,SS)
                W(I,J) = ZERO
                CALL DROT(NCOLS-J+1,W(I-1,J+1),MDW,W(I,J+1),MDW,SC,SS)
  470        CONTINUE
  480     CONTINUE
          MROWS = NCOLS + MVAL + 1

      ELSE
          MROWS = MINPUT
      END IF
C
C      SET THE X(*) ARRAY TO ZERO SO ALL COMPONENTS ARE DEFINED.
      X(1) = ZERO
      CALL DCOPY(NCOLS,X,0,X,1)
C
C     THE ARRAYS IBASIS(*), IBB(*) ARE INITIALIZED BY THE CALLING
C     PROGRAM UNIT.
C     THE COL. SCALING IS DEFINED IN THE CALLING PROGRAM UNIT.
C    'BIG' IS PLUS INFINITY ON THIS MACHINE.
      BIG = D1MACH(2)
      DO 540 J = 1,NCOLS
         ICASE = IND(J)
C     DO CASE(ICASE,4)
         GO TO (490,500,510,520),ICASE

         GO TO 530
C     CASE 1
  490    BU(J) = BIG
         GO TO 530
C     CASE 2
  500    BL(J) = -BIG
         GO TO 530
C     CASE 3
  510    GO TO 530
C     CASE 4
  520    BL(J) = -BIG
         BU(J) = BIG
C     END CASE
  530    CONTINUE
  540 CONTINUE
      DO 560 J = 1,NCOLS
         IF ((BL(J).LE.ZERO.AND.ZERO.LE.BU(J).AND.ABS(BU(J)).LT.
     .       ABS(BL(J))) .OR. BU(J).LT.ZERO) THEN
             T = BU(J)
             BU(J) = -BL(J)
             BL(J) = -T
             SCL(J) = -SCL(J)
             DO 550 I = 1,MROWS
                W(I,J) = -W(I,J)
  550        CONTINUE
         END IF
C
C     INDICES IN SET T(=TIGHT) ARE DENOTED BY NEGATIVE VALUES
C     OF IBASIS(*).
         IF (BL(J).GE.ZERO) THEN
             IBASIS(J) = -IBASIS(J)
             T = -BL(J)
             BU(J) = BU(J) + T
             CALL DAXPY(MROWS,T,W(1,J),1,W(1,NCOLS+1),1)
         END IF

  560 CONTINUE
      NSETB = 0
      ITER = 0
C     END PROCEDURE
      GO TO 50
C     PROCEDURE(PROCESS OPTION ARRAY)
  570 CONTINUE
      IF (IDOPE(5).EQ.1) THEN
          FAC = X(NCOLS+IDOPE(1))
          WT = X(NCOLS+IDOPE(2))
          MVAL = IDOPE(3)

      ELSE
          FAC = 0.
          WT = 1.
          MVAL = 0
      END IF

      ZERO = 0.D0
      ONE = 1.D0
      TWO = 2.D0
      TOLIND =MIN(1.D-5, SQRT(D1MACH(4)))
      TOLSZE = SQRT(D1MACH(4))
C Edit on 950228-1345:
C      ITMAX = 5*MAX(MROWS,NCOLS)
      ITMAX = 5*MAX(MINPUT,NCOLS)
      IPRINT = 0
C
C     CHANGES TO SOME PARAMETERS CAN OCCUR THROUGH THE OPTION
C     ARRAY, IOPT(*).  PROCESS THIS ARRAY LOOKING CAREFULLY
C     FOR INPUT DATA ERRORS.
      LP = 0
      LDS = 0
  580 CONTINUE
      LP = LP + LDS
C
C     TEST FOR NO MORE OPTIONS.
      IP = IOPT(LP+1)
      JP = ABS(IP)
      IF (IP.EQ.99) THEN
          GO TO 590

      ELSE IF (JP.EQ.99) THEN
          LDS = 1
          GO TO 580

      ELSE IF (JP.EQ.1) THEN
C
C     MOVE THE IOPT(*) PROCESSING POINTER.
          IF (IP.GT.0) THEN
              LP = IOPT(LP+2) - 1
              LDS = 0

          ELSE
              LDS = 2
          END IF

          GO TO 580

      ELSE IF (JP.EQ.2) THEN
C
C     CHANGE TOLERANCE FOR RANK DETERMINATION.
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  NERR = 24
                  NCHAR = 89
                  CALL XERRWV(
     .'DBOLSM(). THE OFFSET=(I1) BEYOND POSTION NCOLS=(I2) MUST BE POSIT
     .IVE FOR OPTION NUMBER 2.',NCHAR,NERR,LEVEL,2,IOFF,NCOLS,0,RDUM,
     .                        RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

              TOLIND = X(NCOLS+IOFF)
              IF (TOLIND.LT.D1MACH(4)) THEN
                  NERR = 25
                  NLEVEL = 0
                  NCHAR = 88
                  CALL XERRWV(
     .'DBOLSM(). THE TOLERANCE FOR RANK DETERMINATION=(R1) IS LESS THAN
     .MACHINE PRECISION=(R2).',NCHAR,NERR,NLEVEL,0,IDUM,IDUM,2,TOLIND,
     .                        D1MACH(4))
              END IF

          END IF

          LDS = 2
          GO TO 580

      ELSE IF (JP.EQ.3) THEN
C
C     CHANGE BLOWUP FACTOR FOR ALLOWING VARIABLES TO BECOME
C     INACTIVE.
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  NERR = 26
                  NCHAR = 89
                  CALL XERRWV(
     .'DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST BE POST
     .IVE FOR OPTION NUMBER 3.',NCHAR,NERR,LEVEL,2,IOFF,NCOLS,0,RDUM,
     .                        RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

              TOLSZE = X(NCOLS+IOFF)
              IF (TOLSZE.LE.ZERO) THEN
                  NERR = 27
                  CALL XERRWV(
     .'DBOLSM(). THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING VARI
     .ABLES MUST BE POSITIVE. NOW=(R1).',NCHAR,NERR,LEVEL,0,IDUM,IDUM,1,
     .                        TOLSZE,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

          END IF

          LDS = 2
          GO TO 580

      ELSE IF (JP.EQ.4) THEN
C
C     CHANGE THE MAX. NO. OF ITERATIONS ALLOWED.
          IF (IP.GT.0) THEN
              ITMAX = IOPT(LP+2)
              IF (ITMAX.LE.0) THEN
                  NERR = 28
                  NCHAR = 65
                  CALL XERRWV(
     .'DBOLSM(). THE MAXIMUM NUMBER OF ITERATIONS=(I1) MUST BE POSITIVE.
     .',NCHAR,NERR,LEVEL,1,ITMAX,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

          END IF

          LDS = 2
          GO TO 580

      ELSE IF (JP.EQ.5) THEN
C
C     CHANGE THE FACTOR FOR PRETRIANGULARIZING THE DATA MATRIX.
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              IF (IOFF.LE.0) THEN
                  NERR = 29
                  NCHAR = 89
                  CALL XERRWV(
     .'DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST BE POST
     .IVE FOR OPTION NUMBER 5.',NCHAR,NERR,LEVEL,2,IOFF,NCOLS,0,RDUM,
     .                        RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

              FAC = X(NCOLS+IOFF)
              IF (FAC.LT.ZERO) THEN
                  NERR = 30
                  NLEVEL = 0
                  NCHAR = 104
                  CALL XERRWV(
     .'DBOLSM(). THE FACTOR (NCOLS/MROWS) WHERE PRE-TRIANGULARIZING IS P
     .ERFORMED MUST BE NONNEGATIVE. NOW=(R1).',NCHAR,NERR,NLEVEL,0,IDUM,
     .                        IDUM,1,FAC,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
                  GO TO 600

              END IF

          END IF

          LDS = 2
          GO TO 580

      ELSE IF (JP.EQ.6) THEN
C
C     CHANGE THE WEIGHTING FACTOR (FROM ONE) TO APPLY TO COMPONENTS
C     NUMBERED .GT. MVAL (INITIALLY SET TO 1.)  THIS TRICK IS NEEDED
C     FOR APPLICATIONS OF THIS SUBPROGRAM TO THE HEAVILY WEIGHTED
C     LEAST SQUARES PROBLEM THAT COME FROM EQUALITY CONSTRAINTS.
          IF (IP.GT.0) THEN
              IOFF = IOPT(LP+2)
              MVAL = IOPT(LP+3)
              WT = X(NCOLS+IOFF)
          END IF

          IF (MVAL.LT.0 .OR. MVAL.GT.MINPUT .OR. WT.LE.ZERO) THEN
              NERR = 38
              NLEVEL = 0
              NCHAR = 116
              CALL XERRWV(
     .'DBOLSM(). THE ROW SEPARATOR TO APPLY WEIGHTING (I1) MUST LIE BETW
     .EEN 0 AND MROWS (I2). WEIGHT (R1) MUST BE POSITIVE.',NCHAR,NERR,
     .                    NLEVEL,2,MVAL,MINPUT,1,WT,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
              GO TO 600

          END IF

          LDS = 3
          GO TO 580
C
C     TURN ON DEBUG OUTPUT.
      ELSE IF (JP.EQ.7) THEN
          IF (IP.GT.0) IPRINT = 1
          LDS = 1
          GO TO 580

      ELSE
          NERR = 23
          NCHAR = 46
          CALL XERRWV('DBOLSM. THE OPTION NUMBER=(I1) IS NOT DEFINED.',
     .                NCHAR,NERR,LEVEL,1,IP,IDUM,0,RDUM,RDUM)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 600

      END IF

  590 CONTINUE
C     END PROCEDURE
      GO TO 40
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  600 CONTINUE
      MODE = -NERR
      RETURN
C     END PROCEDURE
C     END PROGRAM
      END
      subroutine dcopy ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DCOPY copies a vector, x, to a vector, y.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    The routine uses unrolled loops for increments equal to one.
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
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Output, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,n

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
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dy(i) = dx(i)
      end do
      if( n .lt. 7 ) return
   40 continue
      do i = m + 1, n, 7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
      end do
      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
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
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in DX.
c
c    Input, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in DY.
c
c    Output, double precision DDOT, the sum of the product of the 
c    corresponding entries of DX and DY.
c
      implicit none

      double precision ddot
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n

      ddot = 0.0d0
      dtemp = 0.0d0
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
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      end do

   60 ddot = dtemp

      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)

c*********************************************************************72

      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
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
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)

c*********************************************************************72

      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
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
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
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
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
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
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
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
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
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
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
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
      SUBROUTINE DMOUT(M,N,LDA,A,IFMT,IDIGIT)
C     REVISED FEB. 27, 1981.
      CHARACTER IFMT*(*)
      CHARACTER*1 ICOL(3)
      DOUBLE PRECISION A(LDA,N)
C
C     DOUBLE PRECISION MATRIX OUTPUT ROUTINE.
C
C  INPUT..
C
C  M,N,LDA,A(*,*) PRINT THE DOUBLE PRECISION ARRAY A(I,J),I=1,...,M,
C                 J=1,...,N, ON OUTPUT UNIT LOUT. LDA IS THE DECLARED
C                 FIRST DIMENSION OF A(*,*) AS SPECIFIED IN THE CALLING
C                 PROGRAM. THE HEADING IN THE FORTRAN FORMAT STATEMENT
C                 IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST STEP.
C                 THE COMPONENTS A(I,J) ARE INDEXED, ON OUTPUT, IN A
C                 PLEASANT FORMAT.
C  IFMT(*)        A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON
C                 OUTPUT UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN
C                 STATEMENT
C                       WRITE(LOUT,IFMT).
C  IDIGIT         PRINT AT LEAST IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
C                 THE SUBPROGRAM WILL CHOOSE THAT INTEGER 6,14,20 OR 28
C                 WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
C                 PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE
C                 UTILIZED TO WRITE EACH LINE OF OUTPUT OF THE ARRAY
C                 A(*,*). (THIS CAN BE USED ON MOST TIME-SHARING
C                 TERMINALS).  IF IDIGIT.GE.0, 133 PRINTING COLUMNS ARE
C                 UTILIZED. (THIS CAN BE USED ON MOST LINE PRINTERS).
C
C  EXAMPLE..
C
C  PRINT AN ARRAY CALLED (SIMPLEX TABLEAU   ) OF SIZE 10 BY 20 SHOWING
C  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
C  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
C
C     DOUBLE PRECISION TABLEU(20,20)
C     M = 10
C     N = 20
C     LDTABL = 20
C     IDIGIT = -6
C     CALL DMOUT(M,N,LDTABL,TABLEU,'(''1SIMPLEX TABLEAU'')',IDIGIT)
C
C
C
C     AUTHORS    JOHN A. WISNIEWSKI   SANDIA LABS ALBUQUERQUE.
C                RICHARD J. HANSON    SANDIA LABS ALBUQUERQUE.
C     DATE       JULY 30,1978.
C
C
C     GET THE UNIT NUMBER WHERE OUTPUTWILL BE WRITTEN.
      J=2
      LOUT=I1MACH(J)
      DATA ICOL(1),ICOL(2),ICOL(3)/'C','O','L'/
      WRITE(LOUT,IFMT)
      IF(M.LE.0.OR.N.LE.0.OR.LDA.LE.0) RETURN
      NDIGIT = IDIGIT
      IF(IDIGIT.EQ.0) NDIGIT = 6
      IF(IDIGIT.GE.0) GO TO 80
C
      NDIGIT = -IDIGIT
      IF(NDIGIT.GT.6) GO TO 20
C
      DO 10 K1=1,N,4
      K2 = MIN0(N,K1+3)
      WRITE(LOUT,1000) (ICOL,I,I=K1,K2)
      DO 10 I=1,M
      WRITE(LOUT,1004) I,(A(I,J),J=K1,K2)
   10 CONTINUE
      RETURN
C
   20 CONTINUE
      IF(NDIGIT.GT.14) GO TO 40
C
      DO 30 K1=1,N,2
      K2 = MIN0(N,K1+1)
      WRITE(LOUT,1001) (ICOL,I,I=K1,K2)
      DO 30 I=1,M
      WRITE(LOUT,1005) I,(A(I,J),J=K1,K2)
   30 CONTINUE
      RETURN
C
   40 CONTINUE
      IF(NDIGIT.GT.20) GO TO 60
C
      DO 50 K1=1,N,2
      K2=MIN0(N,K1+1)
      WRITE(LOUT,1002) (ICOL,I,I=K1,K2)
      DO 50 I=1,M
      WRITE(LOUT,1006) I,(A(I,J),J=K1,K2)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      DO 70 K1=1,N
      K2 = K1
      WRITE(LOUT,1003) (ICOL,I,I=K1,K2)
      DO 70 I=1,M
      WRITE(LOUT,1007) I,(A(I,J),J=K1,K2)
   70 CONTINUE
      RETURN
C
   80 CONTINUE
      IF(NDIGIT.GT.6) GO TO 100
C
      DO 90 K1=1,N,8
      K2 = MIN0(N,K1+7)
      WRITE(LOUT,1000) (ICOL,I,I=K1,K2)
      DO 90 I=1,M
      WRITE(LOUT,1004) I,(A(I,J),J=K1,K2)
   90 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(NDIGIT.GT.14) GO TO 120
C
      DO 110 K1=1,N,5
      K2 = MIN0(N,K1+4)
      WRITE(LOUT,1001) (ICOL,I,I=K1,K2)
      DO 110 I=1,M
      WRITE(LOUT,1005) I,(A(I,J),J=K1,K2)
  110 CONTINUE
      RETURN
C
  120 CONTINUE
      IF(NDIGIT.GT.20) GO TO 140
C
      DO 130 K1=1,N,4
      K2 = MIN0(N,K1+3)
      WRITE(LOUT,1002) (ICOL,I,I=K1,K2)
      DO 130 I=1,M
      WRITE(LOUT,1006) I,(A(I,J),J=K1,K2)
  130 CONTINUE
      RETURN
C
  140 CONTINUE
      DO 150 K1=1,N,3
      K2 = MIN0(N,K1+2)
      WRITE(LOUT,1003) (ICOL,I,I=K1,K2)
      DO 150 I=1,M
      WRITE(LOUT,1007) I,(A(I,J),J=K1,K2)
  150 CONTINUE
      RETURN
 1000 FORMAT(10X,8(5X,3A1,I4,2X))
 1001 FORMAT(10X,5(9X,3A1,I4,6X))
 1002 FORMAT(10X,4(12X,3A1,I4,9X))
 1003 FORMAT(10X,3(16X,3A1,I4,13X))
 1004 FORMAT(1X,'ROW',I4,2X,1P8D14.5)
 1005 FORMAT(1X,'ROW',I4,2X,1P5D22.13)
 1006 FORMAT(1X,'ROW',I4,2X,1P4D28.19)
 1007 FORMAT(1X,'ROW',I4,2X,1P3D36.27)
      END
      function dnrm2 ( n, x, incx )

c*********************************************************************72
c
cc DNRM2 returns the euclidean norm of a vector. 
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c      DNRM2 ( X ) = sqrt ( X' * X )
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
c    Input, double precision X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, double precision DNRM2, the Euclidean norm of X.
c
      implicit none

      integer                           incx, n

      double precision dnrm2 
      double precision                  x( * )


      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      integer               ix
      double precision      absxi, norm, scale, ssq

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
c        call dlassq( n, x, incx, scale, ssq )
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

      dnrm2 = norm

      return
      end
      SUBROUTINE DQEDEV(X,FJ,LDFJ,IGO,IOPT,ROPT)
C***BEGIN PROLOGUE  DQEDEV
C***REFER TO  DQED
C***ROUTINES CALLED  XERROR,CHRCNT
C***END PROLOGUE  DQEDEV
C     REVISED 870204-1100
C     REVISED YYMMDD-HHMM
C     THIS IS THE SUBPROGRAM FOR EVALUATING THE FUNCTIONS
C     AND DERIVATIVES FOR THE NONLINEAR SOLVER, DQED.
C
C     THE USER PROBLEM HAS MCON CONSTRAINT FUNCTIONS,
C     MEQUA LEAST SQUARES EQUATIONS, AND INVOLVES NVARS
C     UNKNOWN VARIABLES.
C
C     WHEN THIS SUBPROGRAM IS ENTERED, THE GENERAL (NEAR)
C     LINEAR CONSTRAINT PARTIAL DERIVATIVES, THE DERIVATIVES
C     FOR THE LEAST SQUARES EQUATIONS, AND THE ASSOCIATED
C     FUNCTION VALUES ARE PLACED INTO THE ARRAY FJ(*,*).
C     ALL PARTIALS AND FUNCTIONS ARE EVALUATED AT THE POINT
C     IN X(*).  THEN THE SUBPROGRAM RETURNS TO THE CALLING
C     PROGRAM UNIT. TYPICALLY ONE COULD DO THE FOLLOWING
C     STEPS:
C
C  IF(IGO.NE.0) THEN
C     STEP 1. PLACE THE PARTIALS OF THE I-TH CONSTRAINT
C             FUNCTION WITH REPECT TO VARIABLE J IN THE
C             ARRAY FJ(I,J), I=1,...,MCON, J=1,...,NVARS.
C  END IF
C     STEP 2. PLACE THE VALUES OF THE I-TH CONSTRAINT
C             EQUATION INTO FJ(I,NVARS+1).
C  IF(IGO.NE.0) THEN
C     STEP 3. PLACE THE PARTIALS OF THE I-TH LEAST SQUARES
C             EQUATION WITH RESPECT TO VARIABLE J IN THE
C             ARRAY FJ(MCON+I,J), I=1,...,MEQUA,
C             J=1,...,NVARS.
C  END IF
C     STEP 4. PLACE THE VALUE OF THE I-TH LEAST SQUARES
C             EQUATION INTO FJ(MCON+I,NVARS+1).
C     STEP 5. RETURN TO THE CALLING PROGRAM UNIT.
C DQEDEV:
C GLOSSARY OF VARIABLES. NOTATION:
C DUMMY-ARG A dummy argument, that is an argument to this prog. unit.
C /S$A$V$E/ SAV Denotes that this variable is local to the routine
C               and is saved between calls to it.
C INTEGER, REAL, DOUBLE PRECISION, LOGICAL, CHARACTER
C               The types of the variables.
C ADJ-ARR An adjustable array, that is an argument to this prog. unit.
C Name      Memory Status  Type     Argument   Uses and comments.
C                                    Status
C ----      -------------  ----     --------   ------------------
C FJ         DUMMY-ARG     REAL      ADJ-ARY
C IGO        DUMMY-ARG     INTEGER
C IOPT       DUMMY-ARG     INTEGER   ADJ-ARY
C LDFJ       DUMMY-ARG     INTEGER
C NERR                     INTEGER
C NMESS                    INTEGER
C ROPT       DUMMY-ARG     REAL      ADJ-ARY
C X          DUMMY-ARG     REAL      ADJ-ARY
C XMESS                    CHAR*128
      DOUBLE PRECISION FJ(LDFJ,*),X(*),ROPT(*)
      INTEGER IGO,IOPT(*)
      CHARACTER XMESS*128
      XMESS =
     . 'DQED. THE EVALUATOR PROGRAM DQEDEV MUST BE WRITTEN BY THE USER.'
      NERR = 09
      IGO = 17
C
C     THE INTENT HERE IS THAT THE EVALUATOR WILL NOT RETURN
C     FROM THE ERROR PROCESSOR CALL.
      CALL CHRCNT(XMESS,NMESS)
      CALL XERROR(XMESS,NMESS,NERR,01)
      RETURN
      END
      SUBROUTINE DQED(DQEDEV,MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,
     .                FNORM,IGO,IOPT,ROPT,IWA,WA)
C***BEGIN PROLOGUE  DQED
C***DATE WRITTEN   851210   (YYMMDD)
C***REVISION DATE  870204   (YYMMDD)
C***CATEGORY NO. K1b,K1b1a2,K1b2a
C***KEYWORDS  NONLINEAR LEAST SQUARES, SIMPLE BOUNDS,
C             LINEAR CONSTRAINTS
C***AUTHOR  HANSON, R. J., SNLA
C           KROGH, F. T., JPL-CIT
C***PURPOSE  SOLVE NONLINEAR LEAST SQUARES AND NONLINEAR
C            EQUATIONS.  USER PROVIDES SIMPLE BOUNDS, LINEAR
C            CONSTRAINTS AND EVALUATION CODE FOR THE FUNCTIONS.
C***ROUTINES CALLED  DQEDMN, DQEDEV, CHRCNT, XERRWV
C***LONG DESCRIPTION
C        SUBROUTINE DQED (DQEDEV, MEQUA, NVARS, MCON, IND, BL, BU, X,
C       *            FJ, LDFJ, RNORM, IGO, IOPT, ROPT,
C       *            IWORK, WORK)
C
C
C  Table of Sections
C  -----------------
C  1. Introduction
C     ------------
C  2. Calling Sequence Explained
C     ------- -------- ---------
C  3. Remarks on the Usage Examples
C     ------- -- --- ----- --------
C  4. Error Messages for DQED()
C     --------------------------
C  5. References
C     ----------
C
C  1. Introduction
C     ------------
C  This  software  package  is  available in both single and double
C  precision.    The   double  precision  version  (type  REAL*8)  is
C  described  below.   For  the  REAL  version  of the
C  documentation  substitute  'REAL' for 'DOUBLE  PRECISION' in the
C  type  statements.  Change the names of the subprograms: 'DQED()'
C  to 'SQED()', 'DQEDEV()' to 'SQEDEV()', and 'D1MACH' to 'R1MACH.'
C
C  The Fortran subprogram, DQED(), solves the constrained nonlinear
C  least squares problem:
C
C  Minimize  the  sum  of  squares  of  MEQUA (generally nonlinear)
C  equations,
C
C       f (x) = 0, I=1,...,MEQUA                    Eq. (1)
C        I
C  where  x  is  a  (vector)  set  of  NVARS unknowns.  (The vector
C  function  with  these  MEQUA  components  is  called f(x) in the
C  discussion  that  follows.)   The components of x may have upper
C  and  lower  bounds  given  by  the  user.   (In  fact all of the
C  possible  cases, no bounds, bounds at one end only, or upper and
C  lower  bounds  can  be  specified.)   Linear  constraints on the
C  unknowns,  more  general than simple bounds,  can also be given.
C  These constraints can be of the equality or inequality type:
C
C       a  x + ... + a       x      =  y , L = 1,...,MCON,
C        L1 1         L,NVARS NVARS     L
C                                                   Eq. (2)
C
C  with  bounds  specified on the y , again given by the user.  The
C                                  L
C  constraints  can  actually  be slightly nonlinear.  In this case
C  the constraints can be described as:
C
C       g (x) =  y , L = 1,...,MCON,                Eq. (2')
C        L        L
C  where bounds are specified on each y .  The functions g (x) must
C                                      L                  L
C  be  defined for all x in the set described by the simple bounds.
C  Experienced users may wish to turn directly to Examples 1 and 2,
C  listed  below,  before  reading  the  subprogram  documentation.
C  There  is  no  size relation required for the problem dimensions
C  MEQUA,  NVARS,  and  MCON  except  that MEQUA and NVARS are both
C  positive,  and   MCON is nonnegative.
C
C  This code package will do a decent job of solving most nonlinear
C  least squares problems that can be expressed as Eqs. (1) and (2)
C  above,  provided  that  continuous  derivatives of the functions
C  with  respect  to the parameters can be computed.  This can also
C  include  problems  where  the derivatives must be computed using
C  some    form    of    numerical    differentiation.    Numerical
C  differentiation is not provided with this software for solving
C  nonlinear  least squares problems.  Refer to the subprogram
C  JACG for numerical differentiation.  (Note: D. Salane has this
C  submitted to TOMS.  It is not included here.)
C
C  The  authors  also  plan  to develop methods that will do a much
C  better  job  of  coping  with  constraints more general than the
C  essentially linear ones indicated above in Eqs. (2)-(2').  There
C  are  nonlinear  least squares problems with innocent looking but
C  highly  nonlinear  constraints  where  this package will fail to
C  work.   The authors also hope to reduce the overhead required by
C  the software.  This high overhead is due primarily to the method
C  used  to  solve  the  inner-loop  quadratic  model problem.  The
C  authors  recommend  that  users consider using the option number
C  14, described below, to suppress use of the quadratic model. The
C  user  may  find  that  the software works quite well without the
C  quadratic  model.  This  may  be important when the function and
C  derivatives  evaluations  are  not expensive but many individual
C  problems are being solved.
C
C  There  are  two  fundamental  ways to use the subprogram DQED().
C  The  most  staightforward way is to make one Fortran CALL to the
C  subprogram  and  obtain  values  for  the unknowns, x.  The user
C  provides  a subprogram DQEDEV(), described below, that gives the
C  subprogram DQED() values of the functions f(x) and g(x), and the
C  derivative  or  Jacobian  matrices  for  f(x)  and  g(x) at each
C  desired  point x.  This usage is called 'forward communication.'
C  An  alternate  way to use the subprogram is to provide an option
C  that  allows  the  user  to communicate these values by 'reverse
C  communication.'   The  subprogram returns to the calling program
C  unit  and  requests  values  for f(x) and g(x), and the Jacobian
C  matrices  for  f(x)  and  g(x)  for  a  given value of x.  (This
C  framework   is   often   required   in  applications  that  have
C  complicated  algorithmic  requirements  for  evaluation  of  the
C  functions.)   An  example  using  both  'forward'  and 'reverse'
C  communication  is  provided  below  (see  Remarks  on  the Usage
C  Examples) for least squares fitting of two exponential functions
C  to five data points.
C
C  2. Calling Sequence Explained
C     ------- -------- ---------
C  There   are  arrays  used  by  the  subprogram  that  must  have
C  dimensions equivalent to the following declarations.
C
C        INTEGER MEQUA, NVARS, MCON, LDFJ, IGO
C        INTEGER IND(NVARS+MCON), IOPT(LIOPT), IWORK(LIWORK)
C
C        DOUBLE PRECISION BL(NVARS+MCON), BU(NVARS+MCON), X(NVARS), RNORM,
C       *ROPT(LROPT), FJ(LDFJ,NVARS+1), WORK(LWORK)
C
C        EXTERNAL DQEDEV
C  The array dimensions must satisfy the bounds:
C        LIOPT .ge.  Number required for options in use.
C        LROPT .ge. Number required for options in use.
C         LDFJ .ge. MEQUA+MCON,
C  The  array  dimensions  for  the arrays IWORK(*) and WORK(*) can
C  change  if either option 14 or option 15 are in use.  For use in
C  the formulas, define:
C
C       MC=MCON
C       ME=MEQUA
C       NV=NVARS
C       MX=MAX(MEQUA,NVARS)
C  If the user is not using option 15, then
C       NT=5.
C  If the user is using option 15, then
C       NT=new number, must be .ge. 2.
C  If the user is not using option 14, then
C       NA=MC+2*NV+NT.
C  If the user is using option 14, then
C       NA=MC+NV+1.
C
C
C  In terms of these values defined above,
C        LIWORK .ge. 3*MC+9*NV+4*NT+NA+10
C         LWORK .ge. NA*(NA+4)+NV*(NT+33)+(ME+MX+14)*NT+9*MC+26
C
C  The  subprogram  DQEDEV  must  be declared in a Fortran EXTERNAL
C  statement:
C
C        EXTERNAL DQEDEV
C
C  Initialize the values of the parameters:
C        MEQUA, NVARS, MCON, IND(*), BL(*), BU(*), X(*), LDFJ,
C        IOPT(*), IWORK(1), IWORK(2),
C
C        CALL DQED  (DQEDEV, MEQUA, NVARS, MCON, IND, BL, BU, X,
C       *            FJ, LDFJ, RNORM, IGO, IOPT, ROPT,
C       *            IWORK, WORK)
C
C  Subprogram parameters:
C
C  DQEDEV (Input)
C  -----
C  This  is  the  name  of  a subprogram that the user will usually
C  supply  for  evaluation  of  the  values  of the constraints and
C  model,  and  the  derivatives  of these functions. The user must
C  provide  this subprogram unless 'reverse communication' is used.
C  A  model  for  writing the subprogram DQEDEV() is provided under
C  the heading Example 1 Using Forward Communication, listed below.
C  Users  may  find  it  convenient to modify this model subprogram
C  when  writing  a  subprogram  for  their  own  application.   If
C  'reverse communication' is used, the user does not need to write
C  a stub or dummy subroutine named DQEDEV().  All that is required
C  is  to  declare exactly this name in an EXTERNAL statement.  The
C  code  package  has a dummy subroutine DQEDEV() that will be used
C  in   the   linking  or  load  step.   Example  2  Using  Reverse
C  Communication, listed below, illustrates this detail.
C
C  MEQUA, NVARS, MCON (Input)
C  ------------------
C  Respectively  they  are:  The number of least squares equations,
C  the  number  of unknowns or variables, and the number of general
C  constraints  for the solution, not including simple bounds.  The
C  values  of  MEQUA  and NVARS must be positive; the value of MCON
C  must  be  nonnegative.   Other  values  for these parameters are
C  errors.
C
C  IND(*),BL(*),BU(*) (Input)
C  ------------------
C  These  arrays  describe  the  form of the simple bounds that the
C  components of x are to satisfy.  Components numbered 1,...,NVARS
C  are  used  to  describe  the  form of the simple bounds that the
C  unknown     are     to     satisfy.      Components     numbered
C  NVARS+1,...,NVARS+MCON  are  used  to  describe  the form of the
C  general  MCON linear constraints.  The first NVARS components of
C  IND(*)  indicate  the type of simple bounds that the solution is
C  to  satisfy.   The  corresponding entries of BL(*) and BU(*) are
C  the bounding value.  The only values of IND(*) allowed are 1,2,3
C  or 4.  Other values are errors.  Specifically:
C
C  IND(J)=1, if x .ge. BL(J) is required; BU(J) is not used.
C                J
C        =2, if x .le. BU(J) is required; BL(J) is not used.
C                J
C        =3, if x .ge. BL(J) and x .le. BU(J) is required.
C                J                J
C        =4, if no bounds on x  are required;
C                             J
C                BL(*),BU(*) are not used.
C  General  linear constraints of the form shown in Eq. (2) require
C  that bounds be given for the linear functions y .  Specifically:
C                                                 L
C
C  IND(NVARS+L)=1,  if y .ge. BL(NVARS+L) is required; BU(*) is not
C                       L
C                 needed.
C
C              =2, if y .le. BU(NVARS+L) is required; BL(*) is not
C                      L
C                  needed.
C              =3, if y .ge. BL(NVARS+L) and y .le. BU(NVARS+L)
C                      L                      L
C
C              =4, if no bounds on y  are required;
C                                   L
C                  BL(*),BU(*) are not used.
C
C  The  values of the bounds for the unknowns may be changed by the
C  user  during  the  evaluation of the functions f(x) and g(x) and
C  their Jacobian matrices.
C
C  X(*),FJ(*,*),LDFJ (Input and Output, except LDFJ which is Input)
C  -----------------
C  The  array  X(*)  contains  the  NVARS  values,  x,   where  the
C  functions  f(x)  and  g(x)  and  their Jacobian matrices will be
C  evaluated  by  the subprogram DQED().  After the computation has
C  successfully  completed, the array X(*) will contain a solution,
C  namely  the  unknowns of the problem, x.  (Success is determined
C  by  an  appropriate  value for IGO.  This parameter is described
C  below.)  Initially  the array X(*) must contain a starting guess
C  for  the  unknowns of the problem, x.  The initial values do not
C  need  to  satisfy the constraints or the bounds.  If they do not
C  satisfy the bounds, then the point will be simply projected onto
C  the  bounds  as  a  first  step.  The first linear change to the
C  values  of  x must satisfy the general constraints.  (It is here
C  that the assumption of their linearity is utilized.)
C
C  The  Fortran  two-dimensional array FJ(*,*) is used to store the
C  linear  constraints  of  Eq. (2) (or more generally the Jacobian
C  matrix  of  the functions g(x) with respect to the variables x),
C  and  the  Jacobian  matrix  of  the function f(x).  The Jacobian
C  matrix of the (linear) constraints is placed in rows 1,...,MCON.
C  The  Jacobian  matrix  of  f(x)  is  placed in rows MCON+1, ...,
C  MCON+MEQUA.   The parameter LDFJ is the leading or row dimension
C  of  the  array  FJ(*,*).  Normally the array FJ(*,*) is assigned
C  values   by   the   user  when  the  nonlinear  solver  requests
C  evaluations  of  the  constraints  g(x)  and  the  function f(x)
C  together  with  the Jacobian matrices G(x) and J(x).  The values
C  of  the  constraint  functions  g (x)  are  placed  in the array
C                                   L
C  FJ(L,NVARS+1),  L=1,...,MCON.  The values of the model functions
C  f (x)  are  placed  in  the array at entries FJ(MCON+I,NVARS+1),
C   I
C  I=1,...,MEQUA.   Note  that the second dimension of FJ(*,*) must
C  be at least NVARS+1 in value.
C
C  RNORM (Output)
C  -----
C  This is the value of the Euclidean length or square root of sums
C  of  squares  of  components  of  the  function  f(x)  after  the
C  approximate solution, x, has been found.  During the computation
C  it  is  updated  and equals the best value of the length of f(x)
C  that has been found.
C
C  IGO (Output; it can be an Input if interrupting the code)
C  ---
C  This  flag  directs  user  action and informs the user about the
C  type  of  results obtained by the subprogram.  The user may find
C  it  convenient  to  treat  the cases abs(IGO) .le. 1 the same as
C  IGO=1.  This has no effect on the solution process.
C
C  The  user  can  interrupt the computation at any time, obtaining
C  the  best  values  of the vector x up to this point,  by setting
C  IGO  to any value .gt. 1 and then return control to DQED().  For
C  example  if  a  calculation  must be done in a certain length of
C  time,  the  user  can,  as  the  end of the time draws near, set
C  IGO=20  and return control to DQED().  It is important that this
C  method  be  used  to  stop further computing, rather than simply
C  proceeding.  The reason for this is that certain flags in DQED()
C  must  be reset before any further solving on subsequent problems
C  can  take  place.  The value of IGO .gt. 1 used to interrupt the
C  computation  is  arbitrary and is the value of IGO returned.  If
C  values of IGO =2,...,18 are used to flag this interrupt, they do
C  not mean the same thing as indicated below.  For this reason the
C  value  IGO=20 is recommended for signaling interrupts in DQED().
C  Another   situation  that  may  occur  is  the  request  for  an
C  evaluation  of  the functions and derivatives at a point x where
C  these can't be evaluated.  If this occurs, set IGO=99 and return
C  control  to  DQED().   This will have the effect of defining the
C  derivatives  to  be  all  zero  and the functions to be 'large.'
C  Thus  a  reduction  in  the trust region around the current best
C  estimate  will occur.  Assigning the value IGO=99 will not cause
C  DQED() to stop computing.
C
C      =0   Place  the  value  of  f(x)  in FJ(MCON+*,NVARS+1).  If
C  'reverse  communication'  is  being used, CALL DQED() again.  If
C  'forward communication' is being used, do a RETURN.
C
C      =1  or  (-1)   Evaluate the Jacobians for the functions g(x)
C  and  f(x) as well as evaluating g(x) and f(x).  Use the vector x
C  that is now in the array X(*) as the values  where this
C  evaluation  will  be performed.  Place the Jacobian matrix
C  for  g(x) in the first MCON rows of FJ(*,*).  Place the Jacobian
C  matrix for f(x) in rows MCON+1,...,MCON+MEQUA in FJ(*,*).  Place
C  the  value of g(x) in FJ(*,NVARS+1).  Place the value of f(x) in
C  FJ(MCON+*,NVARS+1).
C
C  (Users who have complicated functions whose derivatives cannot be
C  computed analytically may want to use the numerical differentiation
C  subroutine JAGC.  This is available on the SLATEC library.)
C
C  If  'reverse communication' is being used, CALL DQED() again.
C  If 'forward communication' is being used, do a RETURN.
C
C  A  value  IGO=(-1)  flags  that  that the number of terms in the
C  quadratic  model  is  being  restricted by the amount of storage
C  given  for  that  purpose.   It  is  suggested,  but  it  is not
C  required,  that  additional  storage  be given for the quadratic
C  model  parameters.   See the following description of The Option
C  Array,  option  number 15, for the way to designate more storage
C  for this purpose.
C
C      =2   The function f(x) has a length less than TOLF.  This is
C  the  value  for  IGO to be expected when an actual zero value of
C  f(x)  is  anticipated.   See the description of The Option Array
C  for the value.
C
C      =3   The  function  f(x)  has  reached a value that may be a
C  local   minimum.   However,  the  bounds  on  the  trust  region
C  defining  the size of the step are being hit at each step.  Thus
C  the  situation  is  suspect.  (Situations of this type can occur
C  when  the  solution  is at infinity in some of the components of
C  the   unknowns,  x.  See the description of The Option Array for
C  ways to avoid this value of output value of IGO.
C
C      =4   The function f(x) has reached a local minimum.  This is
C  the  value  of IGO that is expected when a nonzero value of f(x)
C  is anticipated.  See the description of The Option Array for the
C  conditions that have been satisfied.
C
C      =5   The  model  problem  solver  has  noted a value for the
C  linear or quadratic model problem residual vector length that is
C  .ge.  the  current  value  of  the  function, i.e. the Euclidean
C  length   of  f(x).   This  situation  probably  means  that  the
C  evaluation  of  f(x)  has  more  uncertainty  or  noise  than is
C  possible  to  account for in the tolerances used to note a local
C  minimum.  The value for x is suspect, but a minimum has probably
C  been found.
C
C      =6  A small change (absolute) was noted for the vector x.  A
C  full  model problem step was taken.  The condition for IGO=4 may
C  also  be  satisfied, so that a minimum has been found.  However,
C  this test is made before the test for IGO=4.
C
C      =7   A  small change (relative to the length of x) was noted
C  for  the  vector  x.   A full model problem step was taken.  The
C  condition for IGO=4 may also be satisfied, so that a minimum has
C  been  found.   However,  this  test  is made before the test for
C  IGO=4.
C
C      =8   More  than  ITMAX  iterations  were taken to obtain the
C  solution.   The  value obtained for x is suspect, although it is
C  the   best   set  of  x  values  that  occurred  in  the  entire
C  computation.   See  the  description  of  The  Option  Array for
C  directions  on  how  to  increase  this  value.   (Note that the
C  nominal  value  for ITMAX, 75, is sufficient to solve all of the
C  nonlinear test problems described in Ref. (2).)
C
C      =9-18     Errors  in the usage of the subprogram were noted.
C  The  exact condition will be noted using an error processor that
C  prints  an  informative  message  unless  this printing has been
C  suppressed.   A  minimum  value  has  not been found for x.  The
C  relation  between  IGO  and  the  error number are IGO=NERR + 8.
C  Here  NERR is the identifying number.  See below, Error Messages
C  for DQED().
C  The Option Array
C  --- ------ -----
C  Glossary of Items Modified by Options.  Those items with Nominal
C  Values listed can be changed.
C
C       Names    Nominal         Definition
C                Values
C       -----    -------         ----------
C       FC                       Current value of length of f(x).
C       FB                       Best value of length of f(x).
C       FL                       Value of length of f(x) at the
C                                previous step.
C       PV                       Predicted value of length of f(x),
C                                after the step is taken, using the
C                                approximating model.
C  The  quantity  'eps',  used  below,  is  the  machine  precision
C  parameter.   Its  value  is obtained by a call to the Bell Labs.
C  Port subprogram D1MACH(4).  It is machine dependent.
C                MIN(1.D-5,
C       TOLF     sqrt(eps))      Tolerance for stopping when
C                                FC .le. TOLF.
C                MIN(1.D-5,
C       TOLD     sqrt(eps))      Tolerance for stopping when
C                                change to x values has length
C                                .le. TOLD.
C                MIN(1.D-5,
C       TOLX     sqrt(eps))      Tolerance for stopping when
C                                change to x values has length
C                                .le. TOLX*length of x values.
C       TOLSNR    1.D-5          Tolerance used in stopping
C                                condition IGO=4.  Explained below.
C       TOLP      1.D-5          Tolerance used in stopping
C                                condition IGO=4.  Explained below.
C
C  (The  conditions  (abs(FB-PV).le.TOLSNR*FB  and  abs(FC-PV) .le.
C  TOLP*FB)   and  (ABS(FC-FL).le.TOLSNR*FB) together with  taking
C  a  full  model  step, must be satisfied before the condition  IGO=4
C  is returned.  Decreasing any of the values for  TOLF,  TOLD,  TOLX,
C  TOLSNR,  or  TOLP  will likely increase the number of iterations
C  required for convergence.)
C
C       COND       30.           Largest condition number to allow
C                                when solving for the quadratic
C                                model coefficients.  Increasing
C                                this value may result in more
C                                terms being used in the quadratic
C                                model.
C       TOLUSE   sqrt(eps)       A tolerance that is used to avoid
C                                values of x in the quadratic
C                                model's interpolation of previous
C                                points.  Decreasing this value may
C                                result in more terms being used in
C                                the quadratic model.
C        ITMAX     75            The number of iterations to take
C                                with the algorithm before giving
C                                up and noting it with the value
C                                IGO=8.
C        IPRINT     0            Control the level of printed
C                                output in the solver.  A value
C                                of IPRINT .gt. 0 will result in
C                                output of information about each
C                                iteration.  The output unit used
C                                is obtained using the Bell Labs.
C                                Port subprogram, i. e. I1MACH(2).
C        LEVEL      1            Error processor error level.  See
C                                the SLATEC library documentation
C                                for XERROR() for an explanation.
C        NTERMS     5            One more than the maximum number
C                                of terms used in the quadratic
C                                model.
C
C  IOPT(*) (Input)
C  -------
C  In  order  to  use the option array technique to change selected
C  data within a subprogram, it is necessary to understand how this
C  array  is  processed  within the software.  Let LP designate the
C  processing pointer that moves to positions of the IOPT(*) array.
C  Initially  LP=1,  and  as  each option is noted and changed, the
C  value  of  LP is updated.  The values of IOPT(LP) determine what
C  options get changed.  The amount that LP changes is known by the
C  software  to  be  equal to the value two except for two options.
C  These exceptional cases are the last option (=99) and the 'leap'
C  option  (=13)  which  advances LP by the value in IOPT(LP+1).  A
C  negative  value for IOPT(LP) means that this option is not to be
C  changed.   This aids the programmer in using options;  often the
C  code  for  using  an  option can be in the calling program but a
C  negative value of the option number avoids rewriting code.
C
C  Option Usage Example
C  ------ ----- -------
C  In  the  Fortran code fragment that follows, an example is given
C  where  we  change  the  value  of  TOLF and decrease the maximum
C  number  of  iterations  allowed  from  75  to 30.
C  In this example the dimensions of IOPT(*) and ROPT(*) must
C  satisfy:
C
C        DOUBLE PRECISION ROPT(01)
C        INTEGER IOPT(005)
C        .
C        .
C        .
C  C     SET THE OPTION TO CHANGE THE VALUE OF TOLF.
C
C        IOPT(01)=4
C
C  C     THE NEXT ENTRY POINTS TO THE PLACE IN ROPT(*) WHERE
C  C     THE NEW VALUE OF TOLF IS LOCATED.
C
C        IOPT(02)=1
C  C     THIS IS THE NEW VALUE OF TOLF.  THE SPECIFIC VALUE
C  C     1.D-9 IS USED HERE ONLY FOR ILLUSTRATION.
C
C        ROPT(01)=1.D-9
C
C  C     CHANGE THE NUMBER OF ITERATIONS.
C
C        IOPT(03)=2
C
C  C     THIS NEXT ENTRY IS THE NEW VALUE FOR THE MAXIMUM NUMBER OF
C  C     ITERATIONS.
C
C        IOPT(04)=30
C
C  C     THIS NEXT OPTION IS A SIGNAL THAT THERE ARE NO MORE
C  C     OPTIONS.
C
C        IOPT(05)=99
C        .
C        .
C        .
C        CALL DQED()
C        .
C        .
C        .
C  Option Values   Explanation
C  ------ ------   -----------
C     =99          There are no more options to change.
C                  Normally this is the first and only
C                  option that a user needs to specify,
C                  and it can be simply IOPT(01)=99.  The
C                  total dimension of IOPT(*) must be at
C                  least 17, however.  This can lead to a
C                  hard-to-find program bug if the dimension
C                  is too small.
C
C     = 1          Change the amount of printed output.
C                  The next value of IOPT(*) is the print
C                  level desired, IPRINT.  Any value of
C                  IPRINT .gt. 0 gives all the available
C                  output.
C
C     = 2          Change the value of ITMAX.  The next value
C                  of IOPT(*) is the value of ITMAX desired.
C
C     = 3          Pass prior determined bounds for the box
C                  containing the initial point.  This box is the
C                  trust region for the first move from the initial
C                  point.  The next entry in IOPT(*) points to
C                  the place in ROPT(*) where the NVARS values for
C                  the edges of the box are found.
C
C     = 4          Change the value of TOLF.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLF is found.
C
C     = 5          Change the value of TOLX.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLX is found.
C
C     = 6          Change the value of TOLD.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLD is found.
C
C     = 7          Change the value of TOLSRN.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLSNR is found.
C
C     = 8          Change the value of TOLP.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLP is found.
C
C     = 9          Change the value of TOLUSE.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of TOLUSE is found.
C
C     =10          Change the value of COND.  The next entry of
C                  IOPT(*) points to the place in ROPT(*) where the
C                  new value of COND is found.
C
C     =11          Change the value of LEVEL.  The next entry of
C                  IOPT(*) is the new value of LEVEL.
C
C     =12          Pass an option array to the subprogram DQEDGN()
C                  used as the inner loop solver for the
C                  model problem.  The next entry of IOPT(*) is the
C                  starting location for the option array for
C                  DQEDGN() within the array IOPT(*).  Thus the
C                  option array for DQEDGN() must be a part of
C                  the array IOPT(*).
C
C     =13          Move (or leap) the processing pointer LP for the
C                  option array by the next value in IOPT(*).
C
C     =14          Change a logical flag that suppresses the
C                  use of the quadratic model in the inner
C                  loop.  Use the next value in IOPT(*) for
C                  this flag.  If this value = 1, then never
C                  use the quadratic model.  (Just use the
C                  linear model).  Otherwise, use the quadratic
C                  model when appropriate.  This option decreases
C                  the amount of scratch storage as well as the
C                  computing overhead required by the code package.
C                  A user may want to determine if the application
C                  really requires the use of the quadratic model.
C                  If it does not, then use this option to save
C                  both storage and computing time.
C
C     =15          Change, NTERMS,  the maximum number of array
C                  columns that can be used for saving quadratic
C                  model data.  (The value of NTERMS is one more
C                  than the maximum number of terms used.)  Each
C                  unit increase for NTERMS increases the required
C                  dimension of the array WORK(*) by 2*MEQUA+NVARS.
C                  Use the value in IOPT(LP+1) for the new value
C                  of NTERMS.  Decreasing this value to 2 (its
C                  minimum) decreases the amount of storage
C                  required by the code package.
C
C     =16          Change a logical flag so that 'reverse
C                  communication' is used instead of 'forward
C                  communication.'  Example EX01, listed below,
C                  uses 'forward communication.'  Example EX02,
C                  also listed below, uses 'reverse communication.'
C                  Use the next value in IOPT(*) for
C                  this flag.  If this value = 1, then
C                  use 'reverse communication.'  Otherwise,
C                  use 'forward communication.'  WARNING:  This
C                  usage may not work unless the operating system
C                  saves variables between subroutine calls to DQED.
C
C     =17          Do not allow the flag IGO to return with the
C                  value IGO=3.  This means that convergence will
C                  not be claimed unless a full model step is taken.
C                  Normal output values will then be IGO = 2,4,6 or 7.
C                  Use the next value in IOPT(*) for this flag.  If
C                  this value = 1, then force a full model step.
C                  Otherwise,  do not force a full model step if small
C                  steps are noted.
C
C  IWORK(*), WORK(*) (Input and Output)
C  ----------------
C  These  are  scratch arrays that the software uses for storage of
C  intermediate  results.   It  is  important  not  to  modify  the
C  contents of this storage during the computation.
C
C  The  array  locations  IWORK(1)  and  IWORK(2)  must contain the
C  actual  lengths  of  the  arrays WORK(*) and IWORK(*) before the
C  call to the subprogram.  These array entries are replaced by the
C  actual amount of storage required for each array.  If the amount
C  of  storage  for either array is too small, an informative error
C  message will be printed, and the value IGO=13 or 14 returned.
C
C  The  user may find it useful to let the subprogram DQED() return
C  the  amounts  of storage required for these arrays.  For example
C  set  IWORK(1)=1,  IWORK(2)=1.   The  subprogram will return with
C  IGO=13,     IWORK(1)=required    length    of    WORK(*),    and
C  IWORK(2)=required   length   of  IWORK(*).   (Appropriate  error
C  messages will normally be printed.)
C
C  3. Remarks on the Usage Examples
C     ------- -- --- ----- --------
C  The  following  complete  program units, EX01 and EX02, show how
C  one  can  use  the  nonlinear  solver  for  fitting  exponential
C  functions  to  given data.  These examples are calculations that
C  match  two  terms  of  an  exponential series to five given data
C  points.   There are some subtle points about exponential fitting
C  that   are  important  to  note.     First,  the  signs  of  the
C  exponential arguments are restricted to be nonpositive.
C  The size of the arguments should not be much larger than the start
C  of the time data (reciprocated).  This is the reason the lower
C  bounds are set a bit less than the reciprocal of the time value.
C  In many applications that require exponential modeling this is a
C  natural assumption.  The nonlinear solver allows these bounds
C  on  the arguments explicitly.  In addition, the coefficients are
C  constrained  to  be  nonnegative.   These  bounds  are harder to
C  justify.  The idea is to avoid the situation where a coefficient
C  is  very  large  and negative, and the corresponding exponential
C  argument is also large and negative.  The resulting contribution
C  to  the  series may be very small, but its presence is spurious.
C  Finally,  the  single  general  linear constraint that keeps the
C  arguments  separated  (by  0.05 in this example) is used for two
C  purposes.   First,  it naturally orders these values so that the
C  first  one  is  algebraically  largest.  Second, this constraint
C  moves the parameters from the local minimum corresponding to the
C  initial  values  used  in  the  examples.   This constraint also
C  retains  the  validity of the model function h(t) = w*exp(x*t) +
C  y*exp(z*t).  Namely, if the arguments are allowed to coalesce to
C  the  same value, then the model itself must change.  The form of
C  the model can become h(t)=(a+b*t)*exp(c*t) or h(t) = d*exp(e*t).
C  Either one could occur, and the choice is problem dependent.
C  Example 1  Using Forward Communication
C  ---------  ----- ------- -------------
C      PROGRAM EX01
C
CC     Illustrate the use of the Hanson-Krogh nonlinear least
CC     squares solver for fitting two exponentials to data.
CC
CC     The problem is to find the four variables x(1),...,x(4)
CC     that are in the model function
CC
CC          h(t) = x(1)*exp(x(2)*t) + x(3)*exp(x(4)*t)
CC     There are values of h(t) given at five values of t,
CC     t=0.05, 0.1, 0.4, 0.5, and 1.0.
CC     We also have problem constraints that x(2), x(4) .le. 0, x(1),
CC     x(3) .ge. 0, and a minimal separation of 0.05 between x(2) and
CC     x(4).  Nothing more about the values of the parameters is known
CC     except that x(2),x(4) are approximately .ge. 1/min t.
CC     Thus we have no further knowledge of their values.
CC     For that reason all of the initial values are set to zero.
CC
CC     Dimension for the nonlinear solver.
C      DOUBLE PRECISION FJ(6,5),BL(5),BU(5),X(4),ROPT(001),WA(640)
CC  EDIT on 950228-1300:
C      DOUBLE PRECISION RNORM
C      INTEGER IND(5),IOPT(24),IWA(084)
C
C      EXTERNAL DQEDEX
C
C      DATA LDFJ,LWA,LIWA/6,640,084/
C
C      MCON = 1
C      MEQUA = 5
C      NVARS = 4
CC     Define the constraints for variables.
C      BL(1) = 0.
C      BL(2) = -25.
C      BU(2) = 0.
C      BL(3) = 0.
C      BL(4) = -25.
C      BU(4) = 0.
CC     Define the constraining value (separation) for the arguments.
C      BL(5) = 0.05
CC     Define all of the constraint indicators.
C      IND(1) = 1
C      IND(2) = 3
C      IND(3) = 1
C      IND(4) = 3
C      IND(5) = 1
CC     Define the initial values of the variables.
CC     We don't know anything more, so all variables are set zero.
C      DO 10 J = 1,NVARS
C         X(J) = 0.D0
C   10 CONTINUE
CC     Tell how much storage we gave the solver.
C      IWA(1) = LWA
C      IWA(2) = LIWA
CC     No additional options are in use.
C      IOPT(01) = 99
C      CALL DQED(DQEDEX,MEQUA,NVARS,MCON,IND,BL,BU,X,FJ,LDFJ,RNORM,IGO,
C     .          IOPT,ROPT,IWA,WA)
C      NOUT = 6
C      WRITE (NOUT,9001) (X(J),J=1,NVARS)
C      WRITE (NOUT,9011) RNORM
C      WRITE (NOUT,9021) IGO
C
C      STOP
C
C 9001 FORMAT (' MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))',/,
C     .  ' X(1),X(2),X(3),X(4) = ',/,4F12.6)
C 9011 FORMAT (' RESIDUAL AFTER THE FIT = ',1PD12.4)
C 9021 FORMAT (' OUTPUT FLAG FROM SOLVER =',17X,I6)
C      END
C      SUBROUTINE DQEDEX(X,FJ,LDFJ,IGO,IOPT,ROPT)
CC     This is the subprogram for evaluating the functions
CC     and derivatives for the nonlinear solver, DQED.
CC
CC     The user problem has MCON constraint functions,
CC     MEQUA least squares equations, and involves NVARS
CC     unknown variables.
CC
CC     When this subprogram is entered, the general (near)
CC     linear constraint partial derivatives, the derivatives
CC     for the least squares equations, and the associated
CC     function values are placed into the array FJ(*,*).
CC     All partials and functions are evaluated at the point
CC     in X(*).  Then the subprogram returns to the calling
CC     program unit. Typically one could do the following
CC     steps:
CC
CC     step 1. Place the partials of the i-th constraint
CC             function with respect to variable j in the
CC             array FJ(i,j), i=1,...,MCON, j=1,...,NVARS.
CC     step 2. Place the values of the i-th constraint
CC             equation into FJ(i,NVARS+1).
CC     step 3. Place the partials of the i-th least squares
CC             equation with respect to variable j in the
CC             array FJ(MCON+i,j), i=1,...,MEQUA,
CC             j=1,...,NVARS.
CC     step 4. Place the value of the i-th least squares
CC             equation into FJ(MCON+i,NVARS+1).
CC     step 5. Return to the calling program unit.
C      DOUBLE PRECISION FJ(LDFJ,*),X(*),ROPT(*)
C      DOUBLE PRECISION T(5),F(5)
C      INTEGER IOPT(*)
C
C      DATA T/0.05,0.10,0.40,0.50,1.00/
C      DATA F/2.206D+00,1.994D+00,1.350D+00,1.216D+00,.7358D0/
C
C      DATA MCON,MEQUA,NVARS/1,5,4/
C
CC     Define the derivatives of the constraint with respect to the x(j).
C      FJ(1,1) = 0.D0
C      FJ(1,2) = 1.D0
C      FJ(1,3) = 0.D0
C      FJ(1,4) = -1.D0
CC     Define the value of this constraint.
C      FJ(1,5) = X(2) - X(4)
CC     Define the derivatives and residuals for the data model.
C      DO 10 I = 1,MEQUA
C         E1 = EXP(X(2)*T(I))
C         E2 = EXP(X(4)*T(I))
C         FJ(MCON+I,1) = E1
C         FJ(MCON+I,2) = X(1)*T(I)*E1
C         FJ(MCON+I,3) = E2
C         FJ(MCON+I,4) = X(3)*T(I)*E2
C         FJ(MCON+I,5) = X(1)*E1 + X(3)*E2 - F(I)
C   10 CONTINUE
C      RETURN
C      END
C  Output from Example 1 Program
C  ------ ---- --------- -------
C
C   MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))
C  X(1),X(2),X(3),X(4) =
C      1.999475    -.999801     .500057   -9.953988
C   RESIDUAL AFTER THE FIT =   4.2408D-04
C   OUTPUT FLAG FROM SOLVER =                      4
C
C
C  Example 2  Using Reverse Communication
C  ---------  ----- ------- -------------
C      PROGRAM EX02
C
CC     Illustrate the use of the Hanson-Krogh nonlinear least
CC     squares solver for fitting two exponentials to data.
CC
CC     The problem is to find the four variables x(1),...,x(4)
CC     that are in the model function
CC
CC          h(t) = x(1)*exp(x(2)*t) + x(3)*exp(x(4)*t)
CC     There are values of h(t) given at five values of t,
CC     t=0.05, 0.1, 0.4, 0.5, and 1.0.
CC     We also have problem constraints that x(2), x(4) .le. 0, x(1),
CC     x(3) .ge. 0, and a minimal separation of 0.05 between x(2) and
CC     x(4).  Nothing more about the values of the parameters is known
CC     except that x(2),x(4) are approximately .ge. 1/min t.
CC     Thus we have no further knowledge of their values.
CC     For that reason all of the initial values are set to zero.
CC
CC     Dimension for the nonlinear solver.
C      DOUBLE PRECISION FJ(6,5),BL(5),BU(5),X(4),ROPT(001),WA(640)
CC  EDIT on 950228-1300:
C      DOUBLE PRECISION RNORM
C      INTEGER IND(5),IOPT(24),IWA(084)
C      DOUBLE PRECISION T(5),F(5)
C
C      EXTERNAL DQEDEV
C
C      DATA LDFJ,LWA,LIWA/6,640,084/
C
C      DATA T/0.05,0.10,0.40,0.50,1.00/
C      DATA F/2.206D+00,1.994D+00,1.350D+00,1.216D+00,.7358D0/
C
C      MCON = 1
C      MEQUA = 5
C      NVARS = 4
CC     Define the constraints for variables.
C      BL(1) = 0.
C      BL(2) = -25.
C      BU(2) = 0.
C      BL(3) = 0.
C      BL(4) = -25.
C      BU(4) = 0.
CC     Define the constraining value (separation) for the arguments.
C      BL(5) = 0.05
CC     Define all of the constraint indicators.
C      IND(1) = 1
C      IND(2) = 3
C      IND(3) = 1
C      IND(4) = 3
C      IND(5) = 1
CC     Define the initial values of the variables.
CC     We don't know anything at all, so all variables are set zero.
C      DO 10 J = 1,NVARS
C         X(J) = 0.D0
C   10 CONTINUE
CC     Tell how much storage we gave the solver.
C      IWA(1) = LWA
C      IWA(2) = LIWA
C      NITERS = 0
CC     TELL HOW MUCH STORAGE WE GAVE THE SOLVER.
C      IWA(1) = LWA
C      IWA(2) = LIWA
CC     USE REVERSE COMMUMICATION TO EVALUATE THE DERIVATIVES.
C      IOPT(01)=16
C      IOPT(02)=1
CC     NO MORE OPTIONS.
C      IOPT(03) = 99
C   20 CONTINUE
C      CALL DQED(DQEDEV,MEQUA,NVARS,MCON,IND,BL,BU,X,FJ,LDFJ,RNORM,
C     .IGO,IOPT, ROPT,IWA,WA)
C      IF (IGO.GT.1) GO TO 40
CC     COUNT FUNCTION EVALUATIONS.
C      NITERS = NITERS + 1
CC     DEFINE THE DERIVATIVES OF THE CONSTRAINT WITH RESPECT TO THE X(J).
C      FJ(1,1) = 0.D0
C      FJ(1,2) = 1.D0
C      FJ(1,3) = 0.D0
C      FJ(1,4) = -1.D0
CC     DEFINE THE VALUE OF THIS CONSTRAINT.
C      FJ(1,5) = X(2) - X(4)
CC     DEFINE THE DERIVATIVES AND RESIDUALS FOR THE DATA MODEL.
C      DO 30 I = 1,MEQUA
C          E1 = EXP(X(2)*T(I))
C          E2 = EXP(X(4)*T(I))
C          FJ(MCON+I,1) = E1
C          FJ(MCON+I,2) = X(1)*T(I)*E1
C          FJ(MCON+I,3) = E2
C          FJ(MCON+I,4) = X(3)*T(I)*E2
C          FJ(MCON+I,5) = X(1)*E1 + X(3)*E2 - F(I)
C   30 CONTINUE
C      GO TO 20
C
C   40 CONTINUE
C      NOUT = 6
C      WRITE (NOUT,9001) (X(J),J=1,NVARS)
C      WRITE (NOUT,9011) RNORM
C      WRITE (NOUT,9021) NITERS, IGO
C
C 9001 FORMAT (' MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))',/,
C     . ' X(1),X(2),X(3),X(4) = ',/,4F12.6)
C 9011 FORMAT (' RESIDUAL AFTER THE FIT = ',1PD12.4)
C 9021 FORMAT (' NUMBER OF EVALUATIONS OF PARAMETER MODEL =',I6,/,
C     .          ' OUTPUT FLAG FROM SOLVER =',17X,I6)
C      STOP
C      END
C  Output from Example 2 Program
C  ------ ---- --------- -------
C
C  MODEL IS H(T) = X(1)*EXP(-T*X(2)) + X(3)*EXP(T*X(4))
C  X(1),X(2),X(3),X(4) =
C      1.999475    -.999801     .500057   -9.953988
C   RESIDUAL AFTER THE FIT =   4.2408D-04
C   NUMBER OF EVALUATIONS OF PARAMETER MODEL =    14
C   OUTPUT FLAG FROM SOLVER =                      4
C
C  4. Error Messages for DQED()
C     --------------------------
C   'DQED. VALUE OF MEQUA=NO. OF EQUAS. MUST .GT.0. NOW = (I1).'
C   NERR = 01
C   IGO=9
C
C  'DQED. VALUE OF NVARS=NO. OF EQUAS. MUST .GT.0. NOW = (I1).'
C   NERR = 02
C   IGO=10
C
C  'DQED. VALUE OF MCON=NO. OF EQUAS. MUST .GE.0. NOW = (I1).'
C   NERR = 03
C   IGO=11
C
C  'DQED. INVALID OPTION PROCESSED. I1=IOPT(*) ENTRY.   I2=IOPT(I1).'
C   NERR = 04
C   IGO=12
C
C  'DQED. WA(*) STORAGE SHORT. I1=AMOUNT NEEDED. I2=AMOUNT GIVEN.'
C   NERR = 05
C   IGO=13
C
C  'DQED. IWA(*) STORAGE SHORT. I1=AMOUNT NEEDED. I2=AMOUNT   GIVEN.'
C   NERR = 06
C   IGO=14
C
C  'DQEDMN. INVALID OPTION PROCESSED. I1=IOPT(*) ENTRY.   I2=IOPT(I1).
C
C   NERR=07
C   IGO=15
C
C  'DQEDIP. INVALID OPTION PROCESSED. I1=IOPT(*) ENTRY.  I2=IOPT(I1).'
C
C   NERR=08
C   IGO=16
C
C  'DQED. THE EVALUATOR PROGRAM DQEDEV MUST BE WRITTEN BY THE  USER.'
C   NERR=09
C   IGO=17
C
C  'DQED. BOUND INDICATORS MUST BE 1-4. NOW I1=J, I2=IND(I1).'
C   NERR=10
C   IGO=18
C
C  5. References
C     ----------
C***REFERENCES
C Dongarra, J. J., Bunch, J. R., Moler, C. B, Stewart, G. W.,
C LINPACK User's Guide, Soc. Indust. and Appl. Math, Phil.,
C  PA, (1979).

C Hanson, R. J., "Least Squares with Bounds and Linear
C Constraints," SIAM J. Sci. Stat. Comput., vol. 7, no. 3, July,
C (1986), p. 826-834.

C Schnabel, R. B., Frank, P. D, "Tensor Methods for Nonlinear
C Equations," SIAM J. Num. Anal., vol. 21, no. 5, Oct., (1984),
C p. 815-843.
C***END PROLOGUE  DQED
C     REVISED 870204-1100
C     REVISED YYMMDD-HHMM
      DOUBLE PRECISION BL(*),BU(*),X(*),FJAC(LDFJAC,*)
      DOUBLE PRECISION ROPT(*),WA(*)
      DOUBLE PRECISION FNORM
      REAL RDUM
      INTEGER IND(*),IOPT(*),IWA(*)
      LOGICAL NOQUAD
      CHARACTER XMESS*128
      EXTERNAL DQEDEV
C--PROCEDURES--
C -NAME------TYPE--------ARGS------CLASS-----
C
C  DQEDEV                   6      EXTERNAL
C  DQEDMN                  35      SUBROUTINE
C  CHRCNT                   2      SUBROUTINE
C  XERRWV                  10      SUBROUTINE
C  DQED:
C Name      Memory Status  Type     Argument   Uses and comments.
C                                    Status
C ----      -------------  ----     --------   ------------------
C BL         DUMMY-ARG     REAL      ADJ-ARY Lower bounds
C BU         DUMMY-ARG     REAL      ADJ-ARY Upper bounds
C FJAC       DUMMY-ARG     REAL      ADJ-ARY Jacobian array
C FNORM      DUMMY-ARG     REAL              Norm at solution
C I          /S$A$V$E/ SAV INTEGER           Dummy loop variable
C IDUM       /S$A$V$E/ SAV INTEGER           Dummy for error pack
C IFLAG      /S$A$V$E/ SAV INTEGER           Gate for reverse comm
C IGO        DUMMY-ARG     INTEGER           Directs user action
C IGOOK      /S$A$V$E/ SAV INTEGER           Internal gate for errors
C IIWAAV     /S$A$V$E/ SAV INTEGER           Length claimed for IWA
C IND        DUMMY-ARG     INTEGER   ADJ-ARY Bound indicators
C IOPT       DUMMY-ARG     INTEGER   ADJ-ARY Option array
C IWA        DUMMY-ARG     INTEGER   ADJ-ARY Work array
C IWAAV      /S$A$V$E/ SAV INTEGER           Length claime for WA
C J          /S$A$V$E/ SAV INTEGER           Dummy loop variable
C MILAST     /S$A$V$E/ SAV INTEGER           Last integer in IWA
C MIND       /S$A$V$E/ SAV INTEGER           Point to start IND
C MINDB      /S$A$V$E/ SAV INTEGER           Point to start INDB
C MPJ        /S$A$V$E/ SAV INTEGER           Point to start PJ
C MQC        /S$A$V$E/ SAV INTEGER           Point to start QC
C MUT        /S$A$V$E/ SAV INTEGER           Point to start UT
C MWA        /S$A$V$E/ SAV INTEGER           Point to start WA
C MWJ        /S$A$V$E/ SAV INTEGER           Point to start WJ
C MWLAST     /S$A$V$E/ SAV INTEGER           Last value in WA
C MXB        /S$A$V$E/ SAV INTEGER           Point to start XB
C MXP        /S$A$V$E/ SAV INTEGER           Point to start XP
C MZP        /S$A$V$E/ SAV INTEGER           Point to start ZP
C NALL       /S$A$V$E/ SAV INTEGER           Sum of dimensions
C KP         /S$A$V$E/ SAV INTEGER           Dummy option loop pointer
C LDFJAC     DUMMY-ARG     INTEGER           Row dimension of FJAC
C LEVEL      /S$A$V$E/ SAV INTEGER           Error processor status
C LP         /S$A$V$E/ SAV INTEGER           Dummy option loop pointer
C LPDIFF     /S$A$V$E/ SAV INTEGER           Dummy option loop diff
C MB         /S$A$V$E/ SAV INTEGER           Point to start B
C MBB        /S$A$V$E/ SAV INTEGER           Point to start BB
C MBLB       /S$A$V$E/ SAV INTEGER           Point to start BLB
C MBUB       /S$A$V$E/ SAV INTEGER           Point to start BUB
C MCON       DUMMY-ARG     INTEGER           Number of constraints
C MDX        /S$A$V$E/ SAV INTEGER           Point to start DX
C MDXL       /S$A$V$E/ SAV INTEGER           Point to start DXL
C MEQUA      DUMMY-ARG     INTEGER           Numer of least squares eqs
C MGR        /S$A$V$E/ SAV INTEGER           Point to start GR
C NERR       /S$A$V$E/ SAV INTEGER           Error processor number
C NMESS      /S$A$V$E/ SAV INTEGER           Length of error message
C NOQUAD     /S$A$V$E/ SAV LOGICAL           Flag, suppress quad model
C NPMAX      /S$A$V$E/ SAV INTEGER           Max number of quad terms
C NVARS      DUMMY-ARG     INTEGER           Number of unknowns
C RDUM       /S$A$V$E/ SAV REAL              Dummy variable, error proc
C ROPT       DUMMY-ARG     REAL      ADJ-ARY Option data
C WA         DUMMY-ARG     REAL      ADJ-ARY Working array
C X          DUMMY-ARG     REAL      ADJ-ARY Values of the variables
C XMESS      /S$A$V$E/ SAV CHAR*128          Hold error message
C
      DATA IFLAG/0/
C***FIRST EXECUTABLE STATEMENT  DQED
      ASSIGN 40 TO IGOOK
      IF (IFLAG.EQ.0) THEN
          NOQUAD = .FALSE.
          LEVEL = 1
          IF (MEQUA.LE.0) THEN
              XMESS =
     .      'DQED. VALUE OF MEQUA=NO. OF EQUAS. MUST .GT.0. NOW = (I1).'
              CALL CHRCNT(XMESS,NMESS)
              NERR = 01
              IGO = 9
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,1,MEQUA,IDUM,0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          IF (NVARS.LE.0) THEN
              XMESS =
     .      'DQED. VALUE OF NVARS=NO. OF EQUAS. MUST .GT.0. NOW = (I1).'
              CALL CHRCNT(XMESS,NMESS)
              NERR = 02
              IGO = 10
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,1,NVARS,IDUM,0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          IF (MCON.LT.0) THEN
              XMESS =
     .       'DQED. VALUE OF MCON=NO. OF EQUAS. MUST .GE.0. NOW = (I1).'
              CALL CHRCNT(XMESS,NMESS)
              NERR = 03
              IGO = 11
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,1,MCON,IDUM,0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          DO 10 J = 1,NVARS + MCON
             I = IND(J)
             GO TO (10,10,10,10),I

             XMESS =
     .      'DQED. BOUND INDICATORS MUST BE 1-4. NOW I1=J, I2= IND(I1).'
             CALL CHRCNT(XMESS,NMESS)
             NERR = 10
             IGO = 18
             CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,J,IND(J),0,RDUM,RDUM)
             ASSIGN 50 TO IGOOK

   10     CONTINUE
          NPMAX = 05
C     LOOK THROUGH THE OPTION ARRAY FOR A CHANGE TO NPMAX,
C     THE AMOUNT OF ARRAY STORAGE USED FOR THE QUADRATIC PARAMETERS.
          LP = 1
          LPDIFF = 0
   20     CONTINUE
          LP = LP + LPDIFF
          LPDIFF = 2
          KP = IOPT(LP)
          JP = ABS(KP)
          IF (JP.EQ.99) THEN
              IF (KP.GT.0) GO TO 30
          END IF
C     THIS IS THE ONLY OPTION WHERE THE PROCESSING POINTER
C     MUST BE CHANGED FROM THE VALUE 2.
          IF (JP.EQ.13) LPDIFF = IOPT(LP+1)
C     FOUND A CHANGE TO THE ARRAY SIZE FOR THE QUADRATIC MODEL.
          IF (JP.EQ.15) THEN
              IF (KP.GT.0) NPMAX = IOPT(LP+1)
          END IF
C     SEE IF THE QUADRATIC MODEL IS SUPPRESSED.
C     THIS REQUIRES LESS STORAGE IN THE USER PROGRAM.
          IF (JP.EQ.14) THEN
              IF (KP.GT.0) NOQUAD = IOPT(LP+1) .EQ. 1
          END IF

          IF (JP.LT.1 .OR. JP.GT.17) THEN
C     SAW AN OPTION (OR GARBAGE) THAT IS NOT ON THE LIST.
              XMESS =
     . 'DQED. INVALID OPTION PROOCESSED. I1=IOPT(*) ENTRY. I2=IOPT(I1).'
              NERR = 04
              IGO = 12
              CALL CHRCNT(XMESS,NMESS)
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,LP,IOPT(LP),0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          GO TO 20

   30     CONTINUE
      END IF

      IF (NOQUAD) THEN
          NALL = MCON + NVARS + 2

      ELSE
          NALL = MCON + 2*NVARS + NPMAX + 1
      END IF
C
C     COMPUTE POINTERS INTO WORK SPACE FOR VARIOUS ARRAYS
C     REQUIRED IN MAIN PROGRAMS.
      MDX = 1
      MXB = MDX + NALL + 2
      IF(MCON.GT.0)MXB=MXB+NALL+2
      MB = MXB + NVARS
      MBB = MB + NVARS
      MBLB = MBB + NVARS
      IF(MCON.GT.0)MBLB=MBLB+NALL
      MBUB = MBLB + NALL
      IF(MCON.GT.0)MBUB=MBUB+NALL
      MZP = MBUB + NALL
      MXP = MZP + MEQUA*NPMAX
      MQC = MXP + NVARS*NPMAX
      MWJ = MQC + MAX(MEQUA,NVARS)*NPMAX
      MPJ = MWJ + (NALL)* (NALL+1)
      MGR = MPJ + NVARS + 1
      MDXL = MGR + NVARS
      MWA = MDXL + NVARS + NVARS
      MWLAST = MWA + 9* (MCON+1) + 13* (2*NVARS+NPMAX+1)
C
      MINDB = 3
      MIND = MINDB + NALL + NVARS
      MILAST = MIND + 3* (MCON+1) + 4* (2*NVARS+NPMAX+1)
C     CHECK LENGTHS OF ARRAYS ONCE PER PROBLEM.
      IF (IFLAG.EQ.0) THEN
          IWAAV = IWA(1)
          IIWAAV = IWA(2)
C     RETURN THE ACTUAL AMOUNTS OF STORAGE REQD. FOR WA(*), IWA(*).
          IWA(1) = MWLAST
          IWA(2) = MILAST
          IF (IWAAV.LT.MWLAST) THEN
              XMESS =
     .   'DQED. WA(*) STORAGE SHORT. I1=AMOUNT NEEDED. I2=AMOUNT GIVEN.'
              NERR = 05
              IGO = 13
              CALL CHRCNT(XMESS,NMESS)
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,MWLAST,IWAAV,0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          IF (IIWAAV.LT.MILAST) THEN
              XMESS =
     .  'DQED. IWA(*) STORAGE SHORT. I1=AMOUNT NEEDED. I2=AMOUNT GIVEN.'
              NERR = 06
              IGO = 14
              CALL CHRCNT(XMESS,NMESS)
              CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,MILAST,IIWAAV,0,RDUM,
     .                    RDUM)
              ASSIGN 50 TO IGOOK

          END IF

          IFLAG = 1
      END IF

      GO TO IGOOK

   40 CONTINUE

      CALL DQEDMN(DQEDEV,MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,FNORM,
     .            IGO,IOPT,ROPT,IWA(MIND),WA(MWA),WA(MDX),WA(MXB),
     .            WA(MB),WA(MBB),WA(MBLB),WA(MBUB),IWA(MINDB),NPMAX,
     .            WA(MZP),WA(MXP),WA(MQC),MAX(MEQUA,NVARS),WA(MPJ),
     .            WA(MWJ),NALL,WA(MGR),WA(MDXL))
   50 CONTINUE
      IF (IGO.GT.1) IFLAG = 0
      RETURN
C     TOTAL WORKING STORAGE (FOR DQEDGN) IN WA(*)=
C          9*NALL + 4*NVARS = 9*MCON + 13*NVARS.
C     TOTAL WORKING STORAGE (FOR DQEDGN) IN IWA(*)=
C          3*NALL + NV      = 3*MCON +4*NVARS.
C     IN THE ABOVE FORMULA, NVARS (FOR DQEDGN) IS .LE. 3*NVARS
C     IN TERMS OF THE USER'S VARIABLES.
      END
      SUBROUTINE DQEDGN(MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,FNORM,
     .                  IGO,IOPT,ROPT,IWA,WA)
C***BEGIN PROLOGUE  DQEDGN
C***REFER TO  DQED
C***ROUTINES CALLED  DQEDIP
C***END PROLOGUE  DQEDGN
C  DQEDGN:
C GLOSSARY OF VARIABLES. NOTATION:
C DUMMY-ARG A dummy argument, that is an argument to this prog. unit.
C /S$A$V$E/ SAV Denotes that this variable is local to the routine
C               and is saved between calls to it.
C INTEGER, REAL, DOUBLE PRECISION, LOGICAL, CHARACTER
C               The types of the variables.
C ADJ-ARR An adjustable array, that is an argument to this prog. unit.
C Name      Memory Status  Type     Argument   Uses and comments.
C                                    Status
C ----      -------------  ----     --------   ------------------
C BL         DUMMY-ARG     REAL      ADJ-ARY Model lower bounds
C BU         DUMMY-ARG     REAL      ADJ-ARY Model upper bounds
C FJAC       DUMMY-ARG     REAL      ADJ-ARY Model Jacobian array
C FNORM      DUMMY-ARG     REAL              Model residual norm
C IGO        DUMMY-ARG     INTEGER           direct model action
C IND        DUMMY-ARG     INTEGER   ADJ-ARY Model bound indicators
C IOPT       DUMMY-ARG     INTEGER   ADJ-ARY Option array
C IWA        DUMMY-ARG     INTEGER   ADJ-ARY Working array
C LDFJAC     DUMMY-ARG     INTEGER           Row dim of FJAC(*,*)
C MB         /S$A$V$E/ SAV INTEGER           Pointer to B(*)
C MBB        /S$A$V$E/ SAV INTEGER           Pointer to BB(*)
C MBLB       /S$A$V$E/ SAV INTEGER           Pointer to BLB(*)
C MBUB       /S$A$V$E/ SAV INTEGER           Pointer to BLB(*)
C MCON       DUMMY-ARG     INTEGER           Number, model constraints
C MDX        /S$A$V$E/ SAV INTEGER           Pointer to DX(*)
C MEQUA      DUMMY-ARG     INTEGER           Number, model equations
C MINDB      /S$A$V$E/ SAV INTEGER           Pointer to INDB(*)
C MIWA       /S$A$V$E/ SAV INTEGER           Pointer to IWA(*)
C MWA        /S$A$V$E/ SAV INTEGER           Pointer to WA(*)
C MXB        /S$A$V$E/ SAV INTEGER           Pointer to XB(*)
C NALL       /S$A$V$E/ SAV INTEGER           NVARS+MEQUA
C NVARS      DUMMY-ARG     INTEGER           Number, user variables
C ROPT       DUMMY-ARG     REAL      ADJ-ARY Option array data
C WA         DUMMY-ARG     REAL      ADJ-ARY Working array
C
C     REVISED 870204-1100
C     REVISED YYMMDD-HHMM
      DOUBLE PRECISION BL(*),BU(*),X(*),FJAC(LDFJAC,*),FNORM
      DOUBLE PRECISION ROPT(*),WA(*)
      INTEGER IND(*),IOPT(*),IWA(*)
C     ALLOCATE BLOCKS OF WORKING STORAGE TO LOGICAL ARRAYS.
      NALL = MCON + NVARS
      MDX = 1
      MXB = MDX + 2*NALL + 2
      MB = MXB + NVARS
      MBB = MB + NVARS
      MBLB = MBB + NVARS
      MBUB = MBLB + NALL
      MWA = MBUB + NALL
C
      MINDB = 1
      MIWA = MINDB + NALL + NVARS
      CALL DQEDIP(MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,FNORM,IGO,
     .            IOPT,ROPT,IWA(MIWA),WA(MWA),WA(MDX),WA(MXB),WA(MB),
     .            WA(MBB),WA(MBLB),WA(MBUB),IWA(MINDB))
C     THESE DEFINE THE AMOUNT OF STORAGE FOR THE DOUBLE PRECISION AND
C     INTEGER WORK ARRAYS, WA(*) AND IWA(*).
      MWA = MWA + 6*NVARS + 5*MCON
      MIWA = MIWA + 2*NALL
C     TOTAL WORKING STORAGE IN WA(*)=
C          9*NALL + 4*NVARS = 9*MCON + 13*NVARS.
C     TOTAL WORKING STORAGE IN IWA(*)=
C          3*NALL + NV      = 3*MCON +4*NVARS.
      RETURN
      END
      SUBROUTINE DQEDIP(MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,FB,IGO,
     .                  IOPT,ROPT,IWA,WA,DX,XB,B,BB,BLB,BUB,INDB)
C***BEGIN PROLOGUE  DQEDIP
C***REFER TO  DQED
C***ROUTINES CALLED  DBOCLS,IDAMAX,DCOPY,DNRM2,D1MACH,
C                    CHRCNT,XERRWV
C***END PROLOGUE  DQEDIP
C  DBOCLS                  15      SUBROUTINE
C  IDAMAX    INTEGER        3      FUNCTION
C  DCOPY                    5      SUBROUTINE
C  DNRM2     REAL           3      FUNCTION
C  D1MACH    REAL           1      FUNCTION
C  CHRCNT                   2      SUBROUTINE
C  XERRWV                  10      SUBROUTINE
C DQEDIP:
C GLOSSARY OF VARIABLES. NOTATION:
C DUMMY-ARG A dummy argument, that is an argument to this prog. unit.
C /S$A$V$E/ SAV Denotes that this variable is local to the routine
C               and is saved between calls to it.
C INTEGER, REAL, DOUBLE PRECISION, LOGICAL, CHARACTER
C               The types of the variables.
C ADJ-ARR An adjustable array, that is an argument to this prog. unit.
C Name      Memory Status  Type     Argument   Uses and comments.
C                                    Status
C ----      -------------  ----     --------   ------------------
C ALB        /S$A$V$E/ SAV REAL
C ALFAC      /S$A$V$E/ SAV REAL
C ALPHA      /S$A$V$E/ SAV REAL
C AUB        /S$A$V$E/ SAV REAL
C B          DUMMY-ARG     REAL      ADJ-ARY
C BB         DUMMY-ARG     REAL      ADJ-ARY
C BBOOST     /S$A$V$E/ SAV REAL
C BL         DUMMY-ARG     REAL      ADJ-ARY
C BLB        DUMMY-ARG     REAL      ADJ-ARY
C BOLD       /S$A$V$E/ SAV REAL
C BU         DUMMY-ARG     REAL      ADJ-ARY
C BUB        DUMMY-ARG     REAL      ADJ-ARY
C CHG        /S$A$V$E/ SAV REAL
C CHGFAC     /S$A$V$E/ SAV REAL
C COLNRM     /S$A$V$E/ SAV REAL
C C1516      /S$A$V$E/ SAV REAL
C DX         DUMMY-ARG     REAL      ADJ-ARY
C IPLS       /S$A$V$E/ SAV INTEGER
C IPRINT     /S$A$V$E/ SAV INTEGER
C ITERS      /S$A$V$E/ SAV INTEGER
C ITMAX      /S$A$V$E/ SAV INTEGER
C IWA        DUMMY-ARG     INTEGER   ADJ-ARY
C J          /S$A$V$E/ SAV INTEGER
C JP         /S$A$V$E/ SAV INTEGER
C K          /S$A$V$E/ SAV INTEGER
C KL         /S$A$V$E/ SAV INTEGER
C KP         /S$A$V$E/ SAV INTEGER
C LDFJAC     DUMMY-ARG     INTEGER
C LEVEL      /S$A$V$E/ SAV INTEGER
C  EDIT on 950228-1300. REMOVE references to:
C LINPRB     /S$A$V$E/ SAV INTEGER
C LP         /S$A$V$E/ SAV INTEGER
C LPDIFF     /S$A$V$E/ SAV INTEGER
C MCON       DUMMY-ARG     INTEGER
C MEQUA      DUMMY-ARG     INTEGER
C MODE       /S$A$V$E/ SAV INTEGER
C NALL       /S$A$V$E/ SAV INTEGER
C NERR       /S$A$V$E/ SAV INTEGER
C NEWBST     /S$A$V$E/ SAV LOGICAL
C NEWOPT     /S$A$V$E/ SAV LOGICAL
C NMESS      /S$A$V$E/ SAV INTEGER
C NOUT       /S$A$V$E/ SAV INTEGER
C NVARS      DUMMY-ARG     INTEGER
C ONE        /S$A$V$E/ SAV REAL
C PASSB      /S$A$V$E/ SAV LOGICAL
C DXNRM      /S$A$V$E/ SAV REAL
C FB         DUMMY-ARG     REAL
C FC         /S$A$V$E/ SAV REAL
C FJAC       DUMMY-ARG     REAL      ADJ-ARY
C FL         /S$A$V$E/ SAV REAL
C FULNWT     /S$A$V$E/ SAV LOGICAL
C GVAL       /S$A$V$E/ SAV REAL
C ICASE      /S$A$V$E/ SAV INTEGER
C IFLAG      /S$A$V$E/ SAV INTEGER
C IGO        DUMMY-ARG     INTEGER
C IGOIOV     /S$A$V$E/ SAV INTEGER
C IGOPOA     /S$A$V$E/ SAV INTEGER
C IGOTFC     /S$A$V$E/ SAV INTEGER
C IGOTNC     /S$A$V$E/ SAV INTEGER
C IND        DUMMY-ARG     INTEGER   ADJ-ARY
C INDB       DUMMY-ARG     INTEGER   ADJ-ARY
C IOPT       DUMMY-ARG     INTEGER   ADJ-ARY
C PB         /S$A$V$E/ SAV REAL
C PD         /S$A$V$E/ SAV REAL
C PV         /S$A$V$E/ SAV REAL
C RB         /S$A$V$E/ SAV REAL
C RDUM       /S$A$V$E/ SAV REAL
C RETREA     /S$A$V$E/ SAV LOGICAL
C RG         /S$A$V$E/ SAV REAL
C RNORMC     /S$A$V$E/ SAV REAL
C ROPT       DUMMY-ARG     REAL      ADJ-ARY
C SEMIBG     /S$A$V$E/ SAV REAL
C T          /S$A$V$E/ SAV REAL
C TERM       /S$A$V$E/ SAV LOGICAL
C TOLD       /S$A$V$E/ SAV REAL
C TOLF       /S$A$V$E/ SAV REAL
C TOLP       /S$A$V$E/ SAV REAL
C TOLSNR     /S$A$V$E/ SAV REAL
C TOLX       /S$A$V$E/ SAV REAL
C TWO        /S$A$V$E/ SAV REAL
C T2         /S$A$V$E/ SAV REAL
C WA         DUMMY-ARG     REAL      ADJ-ARY
C X          DUMMY-ARG     REAL      ADJ-ARY
C XB         DUMMY-ARG     REAL      ADJ-ARY
C XMESS      /S$A$V$E/ SAV CHAR*128
C ZERO       /S$A$V$E/ SAV REAL
C
C     REVISED 870204-1100
C     REVISED YYMMDD-HHMM
      DOUBLE PRECISION BL(*),BU(*),X(*),FJAC(LDFJAC,*)
      DOUBLE PRECISION DX(*),XB(*),B(*),BB(*)
      DOUBLE PRECISION BLB(*),BUB(*),ROPT(*),WA(*)
      DOUBLE PRECISION ALB,ALFAC,ALPHA,AUB,BBOOST
      DOUBLE PRECISION BOLD,CHG,CHGFAC,COLNRM
      DOUBLE PRECISION DXNRM,C1516,FB,FC,FL
      DOUBLE PRECISION RG,RB,PB,PD
      DOUBLE PRECISION GVAL,ONE,PV
      DOUBLE PRECISION RNORMC,SEMIBG,T,TOLD,TOLF,TOLP
      DOUBLE PRECISION TOLSNR,TOLX,TWO,T2,ZERO
      DOUBLE PRECISION D1MACH,DNRM2
      REAL RDUM
      INTEGER IND(*),INDB(*),IOPT(*),IWA(*)
      CHARACTER XMESS*128
C     OPTIONS..
C     1    SET THE PRINTED OUTPUT OFF/ON.  REQUIRES TWO ENTRIES
C          IN IOPT(*).  IF IOPT(*+1)=0, NO PRINTING; =1 PRINT.
C          (THIS PRINTOUT SHOWS VARIOUS QUANTITIES ASSOCIATED
C          WITH THE NONLINEAR PROBLEM BEING SOLVED. GOOD FOR
C          DEBUGGING A PROBLEM IN CASE OF TROUBLE.
C     2    SET THE MAXIMUM NUMBER OF INTERATIONS.  REQUIRES TWO ENTRIES
C          IN IOPT(*).  USE IOPT(*+1)= MAXIMUM NUMBER OF ITERATIONS.
C     3    PASS INITIAL BOUNDS FOR THE TRUST REGION TO THE NONLINEAR
C          SOLVER. REQUIRES TWO ENTRIES IN IOPT(*). USE IOPT(*+1) AS A
C          POINTER INTO ROPT(*) FOR START OF THE NVARS BOUNDS.
C  EDIT on 950228-1300:
C      LOGICAL RETREA,TERM,FULNWT,PASSB,NEWBST,NEWOPT,LINPRB
       LOGICAL RETREA,TERM,FULNWT,PASSB,NEWBST,NEWOPT

      DATA IFLAG/0/
C--PROCEDURES--
C -NAME------TYPE--------ARGS------CLASS-----
C
C  DBOCLS                  15      SUBROUTINE
C  IDAMAX    INTEGER        3      FUNCTION
C  DCOPY                    5      SUBROUTINE
C  DNRM2     REAL           3      FUNCTION
C  D1MACH    REAL           1      FUNCTION
C  CHRCNT                   2      SUBROUTINE
C  XERRWV                  10      SUBROUTINE
      GO TO (50),IFLAG

      ZERO = 0.D0
      ONE = 1.D0
      TWO = 2.D0
C     DO(PROCESS OPTION ARRAY)
      ASSIGN 10 TO IGOPOA
      GO TO 470

   10 CONTINUE
C     DO(INITIALIZE OTHER VALUES)
      ASSIGN 20 TO IGOIOV
      GO TO 450

   20 CONTINUE
C     SET SO X(*)-DX(*) UPDATE WILL BE CORRECT FIRST TIME.
      DX(1) = ZERO
      CALL DCOPY(NVARS,DX,0,DX,1)
      K = 0
C     D1MACH(2)="INFINITY" ON THIS MACHINE.
      FB = D1MACH(2)
      DXNRM = FB
      FL = ZERO
C
C     LINEAR PROBLEM RESIDUAL.
      PV = ZERO
      RETREA = .FALSE.
      FULNWT = .FALSE.
      TERM = .FALSE.
   30 CONTINUE
      IF ( .NOT. RETREA) ITERS = ITERS + 1
      IF (RETREA) THEN
C     MUST RETREAT TO BEST X VALUES.
          K = 0
          KL = -1
          FL = FB
          CALL DCOPY(NVARS,XB,1,X,1)

      ELSE
          KL = K
          DO 40 J = 1,NVARS
             X(J) = X(J) - DX(J)
   40     CONTINUE
          IF (TERM) THEN
              IFLAG = 0
              GO TO 390

          END IF

      END IF

      IGO = 1
      IFLAG = 1
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 440

   50 CONTINUE
      FC = DNRM2(MEQUA,FJAC(MCON+1,NVARS+1),1)
C     DO(TEST FOR CONVERGENCE)
      ASSIGN 60 TO IGOTFC
      GO TO 400

   60 CONTINUE
      IF (TERM) THEN
          IFLAG = 0
          GO TO 390

      END IF

      NEWBST = FC .LT. FB .OR. (MCON.GT.0 .AND.ITERS.EQ.2)
      IF (NEWBST) K = 0
      IF (K.EQ.0) THEN
          RG = ZERO
          PB = ZERO
          PD = D1MACH(2)
C     WANT TO POSITION AT BEST X VALUES.
          FB = FC
C     DO CASE(2-KL,3)
          GO TO (70,90,110),2 - KL

          GO TO 120

   70     CONTINUE
C     CASE 1
C     IMMEDIATELY GOT A NEW BEST X.
          IF (T2.LE.0.25D0) THEN
              BBOOST = ONE
              CHG = MAX(4.D0*T2,.1D0)
          END IF

          DO 80 J = 1,NVARS
             BB(J) = CHG*BB(J)
   80     CONTINUE
C     THIS CODE FOR ALPHA HAS THE FOLLOWING EFFECT.
C     IF FC .EQ. PV, THEN ALPHA=ALFAC.
C     IF FC**2 .EQ. PV*FL THEN ALPHA=2.-1./ALFAC
C     IF FC**2 IS MUCH LARGER THAN PV*FL, THEN ALPHA=1.
          T = FC - PV
          IF (T.EQ.ZERO) THEN
              ALPHA = ALFAC

          ELSE
              ALPHA = (PV* (FL-PV))/ (FC+PV)/ (ALFAC-ONE)
              ALPHA = (ABS(T)+ALFAC*ALPHA)/ (ABS(T)+ALPHA)
          END IF

          ALFAC = 1.5D0*ALPHA
          BBOOST = MIN(1.5D0*ALPHA*BBOOST,SEMIBG)
          GO TO 140

   90     CONTINUE
C     CASE 2
C     AT THE INITIAL X.
          ALFAC = 256.D0
          DO 100 J = 1,NVARS
             IF ( .NOT. PASSB) BB(J) = -X(J)
             IF (BB(J).EQ.ZERO) THEN
                 COLNRM = DNRM2(MEQUA,FJAC(MCON+1,J),1)
                 IF (COLNRM.NE.ZERO) BB(J) = -FC/COLNRM
             END IF

             IF (BB(J).EQ.ZERO) BB(J) = -ONE
             XB(J) = X(J)
             B(J) = BB(J)
  100     CONTINUE
          ALPHA = ONE
          BBOOST = 0.5D0
C     EXIT IF
          GO TO 170

  110     CONTINUE
C     CASE 3
C     RETREAT TO BEST X.
          IF (ALFAC.NE.256.D0) THEN
              ALPHA = MIN(ONE/ALFAC,.25D0)
              ALFAC = 1.25D0

          ELSE
              ALPHA = .25D0*ALPHA
          END IF

          BBOOST = .25D0
          GO TO 140

  120     CONTINUE
C     CASE OTHER
C     NOT IMMEDIATELY A BEST X.
          RB = ZERO
          DO 130 J = 1,NVARS
             RB = MAX(RB,ABS((XB(J)-X(J))/BB(J)))
  130     CONTINUE
          ALPHA = RB
          ALFAC = TWO
          BBOOST = (8.D0/7.D0+RG)/ (2.D0/7.D0+RG)
C     END CASE
  140     CONTINUE
          DO 150 J = 1,NVARS
             DX(J) = XB(J) - X(J)
             IF (DX(J).EQ.ZERO) THEN
                 B(J) = ALPHA*BB(J)

             ELSE
                 XB(J) = X(J)
                 B(J) = SIGN(ALPHA*BB(J),DX(J)) + BBOOST*DX(J)
             END IF

            BB(J) = SIGN(MIN(SQRT(D1MACH(2)),ABS(B(J))),B(J))
  150     CONTINUE

      ELSE
C     COMPUTE A GAUGE FOR RETREATING IF DID NOT
C     GET A NEW BEST.
          IF (K.EQ.1) THEN
              PB = PV
              PD = 1.5D0* (FB+PB* (PB/FB)) - 4.D0*PB
          END IF

          ALPHA = (.5D0*FC+FL)/ (FC+FL)
          CHG = MIN(ALPHA*CHG,T2)
          CHG=MAX(CHG,.1D0)
          DO 160 J = 1,NVARS
             B(J) = CHG*B(J)
             IF (K.EQ.1) BB(J) = B(J)
  160     CONTINUE
      END IF

  170 CONTINUE
C     DO(TEST FOR CONVERGENCE)
      ASSIGN 180 TO IGOTFC
      GO TO 400

  180 CONTINUE
      IF (TERM) THEN
          IFLAG = 0
          GO TO 390

      END IF

      K = K + 1
C     SOLVE LINEAR BOUNDED PROBLEM.
      DO 240 J = 1,NVARS
         IF (B(J).LT.ZERO) THEN
             ALB = B(J)
             IF (DX(J).EQ.ZERO) THEN
C     THIS CASE IS REQD. TO AVOID USING BUB(*) AT THE INITIAL PT.
                 AUB = -C1516*ALB

             ELSE
                 AUB = MIN(-C1516*ALB,-DX(J)+BUB(J))
             END IF

         ELSE
             AUB = B(J)
             IF (DX(J).EQ.ZERO) THEN
                 ALB = -C1516*AUB

             ELSE
                 ALB = MAX(-C1516*AUB,-DX(J)+BLB(J))
             END IF

         END IF
C Edit on 950228-1300:
C      IF(LINPRB) THEN
C          AUB=D1MACH(2)
C          ALB=-AUB
C      END IF

C     RESTRICT THE STEP FURTHER IF USER GIVES BOUNDS.
         ICASE = IND(J)
C     DO CASE(ICASE,4)
         GO TO (190,200,210,220),ICASE

  190    CONTINUE
C     CASE 1
         AUB = MIN(AUB,X(J)-BL(J))
         GO TO 230

  200    CONTINUE
C     CASE 2
         ALB = MAX(ALB,X(J)-BU(J))
         GO TO 230

  210    CONTINUE
C     CASE 3
         AUB = MIN(AUB,X(J)-BL(J))
         ALB = MAX(ALB,X(J)-BU(J))
         GO TO 230

  220    CONTINUE
C     CASE 4
C     END CASE
  230    CONTINUE
         BLB(J) = ALB
C     THIS NEXT LINE IS TO GUARANTEE THAT THE LOWER BOUND
C     IS .LE. THE UPPER BOUND.
         AUB = MAX(AUB,ALB)
         BUB(J) = AUB
         INDB(J) = 3
  240 CONTINUE
C     SEE IF USER HAS GIVEN GENERAL CONSTRAINTS.
      DO 300 J = NVARS + 1,NALL
         ICASE = IND(J)
         GVAL = FJAC(J-NVARS,NVARS+1)
C     DO CASE(ICASE,4)
         GO TO (250,260,270,280),ICASE

  250    CONTINUE
C     CASE 1
         BLB(J) = - (GVAL-BL(J))
         INDB(J) = 1
         GO TO 290

  260    CONTINUE
C     CASE 2
         BUB(J) = - (GVAL-BU(J))
         INDB(J) = 2
         GO TO 290

  270    CONTINUE
C     CASE 3
         BLB(J) = - (GVAL-BL(J))
         BUB(J) = - (GVAL-BU(J))
         INDB(J) = 3
         GO TO 290

  280    CONTINUE
C     CASE 4
         INDB(J) = 4
C     END CASE
  290    CONTINUE
  300 CONTINUE
C     SOLVE THE LEAST SQUARES PROBLEM WITH BOUNDS AND LINEAR
C     CONSTRAINTS.  THESE BOUNDS CAN COME FROM THE USER OR
C     THE ALGORITHM.
      CALL DBOCLS(FJAC,LDFJAC,MCON,MEQUA,NVARS,BLB,BUB,INDB,IOPT(IPLS),
     .            DX,RNORMC,PV,MODE,WA,IWA)
      IF (IPRINT.GT.0) THEN
          WRITE (NOUT,9011) ITERS,FC,PV,K,KL,FB,ALPHA,BBOOST
          WRITE (NOUT,9001) '  +X=', (X(J),J=1,NVARS)
          WRITE (NOUT,9001) ' +XB=', (XB(J),J=1,NVARS)
          WRITE (NOUT,9001) ' +DX=', (DX(J),J=1,NALL)
          WRITE (NOUT,9001) ' + B=', (B(J),J=1,NVARS)
          WRITE (NOUT,9001) ' +LB=', (BLB(J),J=1,NALL)
          WRITE (NOUT,9001) ' +UB=', (BUB(J),J=1,NALL)
          WRITE (NOUT,'('' +///END OF ITERATION.///'')')
      END IF
C     TEST FOR NOISE IN LINEAR PROBLEM SOLN.
      TERM=MCON.EQ.0.AND.(PV.GE.FC)
      TERM=.FALSE.
      IF (TERM) THEN
          IF (IPRINT.GT.0) THEN
              WRITE (NOUT,9021) PV,FC
          END IF

          CALL DCOPY(NVARS,XB,1,X,1)
          IGO = 4
          IFLAG = 0
          GO TO 390

      END IF

      RG = MAX(RG, (PV-PB)/PD)
      IF ( .NOT. RETREA) THEN
          CHG = ONE
          T2 = ZERO
          DO 360 J = 1,NVARS
             BOLD = B(J)
             T = DX(J)/BOLD
C    IF USER GIVES BOUNDS, AND THESE BOUNDS ARE HIT,
C    DO NOT DETRACT FROM DECLARING A FULL NEWTON STEP.
             ICASE = IND(J)
             GO TO (310,320,330,340),ICASE
C     CASE 1
  310        ALB = (X(J)-BL(J))/BOLD
             AUB = -SEMIBG
             GO TO 350
C     CASE 2
  320        AUB = (X(J)-BU(J))/BOLD
             ALB = -SEMIBG
             GO TO 350
C     CASE 3
  330        ALB = (X(J)-BL(J))/BOLD
             AUB = (X(J)-BU(J))/BOLD
             GO TO 350
C     CASE 4
  340        CONTINUE
             ALB = -SEMIBG
             AUB = -SEMIBG
C     END CASE
  350        CONTINUE
             IF (T.EQ.ONE) THEN
                 T2 = ONE
                 B(J) = BOLD + BOLD
                 CHG = CHG*CHGFAC

             ELSE
                 IF (ABS(T).LT..25D0 .AND. DX(J).NE.ZERO) THEN
                     B(J) = SIGN(.25D0*BOLD,DX(J)) + 3.D0*DX(J)

                 ELSE
                     B(J) = SIGN(BOLD,DX(J))
                 END IF

             END IF
C     THIS TEST AVOIDS THE USER BOUNDS IN DECLARING A NEWTON STEP.
             IF (ABS(ALB-T).GE..01D0*ABS(T) .AND. ABS(AUB-T).GE..01D0*
     .           ABS(T)) THEN
                 IF (T.GT.ZERO) THEN
                     T2 = MAX(T2,T)

                 ELSE
                     T2 = MAX(T2,-T/C1516)
                 END IF

             END IF

  360     CONTINUE
          FULNWT = T2 .LT. .99D0
          FL = FC
          DXNRM = ABS(DX(IDAMAX(NVARS,DX,1)))
C     DO BLOCK
C     TEST FOR SMALL ABSOLUTE CHANGE IN X VALUES.
          TERM = DXNRM .LT. TOLD .AND. FULNWT
          IF (TERM) THEN
              IGO = 5
C     EXIT BLOCK
              GO TO 370

          END IF

          TERM = DXNRM .LT. DNRM2(NVARS,X,1)*TOLX .AND. FULNWT
          TERM = TERM .AND. (ITERS.GT.1)
          IF (TERM) THEN
              IGO = 6
C     EXIT BLOCK
              GO TO 370

          END IF
C     EXIT IF
          GO TO 380
C     END BLOCK
  370     CONTINUE
          GO TO 30

      END IF

  380 CONTINUE
      GO TO 30

  390 CONTINUE
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 440
C     PROCEDURE(TEST FOR CONVERGENCE)
  400 CONTINUE
C     TEST FOR SMALL FUNCTION NORM.
      TERM = FC .LE. TOLF .OR. TERM
C     IF HAVE CONSTRAINTS MUST ALLOW AT LEAST ONE MOVE.
      TERM = TERM .AND. (MCON.EQ.0 .OR. ITERS.GT.1)
      IF (TERM) THEN
          IGO = 2
C     EXIT PROCEDURE
          GO TO 420

      END IF
C     DO(TEST FOR NO CHANGE)
      ASSIGN 410 TO IGOTNC
      GO TO 430

  410 CONTINUE
      TERM = TERM .AND. .NOT. RETREA
      IF (TERM) THEN
          IGO = 3
C     EXIT PROCEDURE
          GO TO 420

      END IF

      TERM = ITERS .GE. ITMAX
      IF (TERM) THEN
          IGO = 7
      END IF
C     END PROCEDURE
  420 CONTINUE
      GO TO IGOTFC
C     PROCEDURE(TEST FOR NO CHANGE)
  430 CONTINUE
      T = SQRT(MAX(ZERO, (FL-PV)* (FL+PV)))
      TERM = (ABS(FB-FC).LE.TOLSNR*FB) .AND. (T.LE.PV*TOLP)
      TERM = TERM .AND. (ABS(FC-FL).LE.FB*TOLSNR)
      TERM = TERM .AND. FULNWT
C     END PROCEDURE
      GO TO IGOTNC
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
  440 CONTINUE
      RETURN
C     END PROCEDURE
C     PROCEDURE(INITIALIZE OTHER VALUES)
  450 CONTINUE
C Edit on 950228-1300:
C      LINPRB = IGO.EQ.0
      ITERS = 0
      NALL = MCON + NVARS
      CHGFAC = TWO** (-ONE/REAL(NVARS))
      C1516 = 15.D0/16.D0
      SEMIBG = 1.D10
C     MAKE SURE THAT VARIABLES SATISFY THE BOUNDS AND CONSTRAINTS.
      DO 460 J = 1,NALL
         BLB(J) = BL(J)
         BUB(J) = BU(J)
         INDB(J) = IND(J)
  460 CONTINUE
C     END PROCEDURE
      GO TO IGOIOV
C     PROCEDURE(PROCESS OPTION ARRAY)
  470 CONTINUE
      IPRINT = 0
C     I1MACH(2)=STANDARD OUTPUT UNIT.
      NOUT = I1MACH(2)
C     D1MACH(4)=RELPR=MACHINE REL. PREC.
      T = D1MACH(4)
      TOLF = T
      TOLX = TOLF
      TOLD = TOLF
      TOLSNR = 1.D-3
      TOLP = 1.D-3
      ITMAX = 18
      PASSB = .FALSE.
      LEVEL = 1
      IPLS = 0
      LPDIFF = 0
      LP = 1
  480 CONTINUE
      LP = LP + LPDIFF
      LPDIFF = 2
      KP = IOPT(LP)
      NEWOPT = KP .GT. 0
      JP = ABS(KP)
C     SEE IF THIS IS THE LAST OPTION..
      IF (JP.EQ.99) THEN
          IF (NEWOPT) THEN
C     THE POINTER TO THE START OF OPTIONS FOR THE LINEAR
C     SOLVER MUST SATISFY THE REQUIREMENTS FOR THAT OPTION ARRAY.
              IF (IPLS.EQ.0) IPLS = LP
              GO TO 490

          ELSE
              LPDIFF = 1
              GO TO 480

          END IF

      END IF
C     CHANGE PRINT OPTION.
      IF (JP.EQ.1) THEN
          IF (NEWOPT) IPRINT = IOPT(LP+1)
          GO TO 480

      END IF
C     SEE IF MAX. NUMBER OF ITERATIONS CHANGING.
      IF (JP.EQ.2) THEN
          IF (NEWOPT) ITMAX = IOPT(LP+1)
          GO TO 480

      END IF
C     SEE IF BOUNDS FOR THE TRUST REGION ARE BEING PASSED.
      IF (JP.EQ.3) THEN
          IF (NEWOPT) THEN
              CALL DCOPY(NVARS,ROPT(IOPT(LP+1)),1,BB,1)
              PASSB = .TRUE.
          END IF

          GO TO 480

      END IF
C     CHANGE TOLERANCE ON THE LENGTH OF THE RESIDUALS.
      IF (JP.EQ.4) THEN
          IF (NEWOPT) TOLF = ROPT(IOPT(LP+1))
          GO TO 480

      END IF
C     CHANGE TOLERANCE ON THE NORM OF THE RELATIVE
C     CHANGE TO THE PARAMETERS.
      IF (JP.EQ.5) THEN
          IF (NEWOPT) TOLX = ROPT(IOPT(LP+1))
          GO TO 480

      END IF
C     CHANGE TOLERANCE ON ABSOLUTE CHANGE TO THE PARAMETERS.
      IF (JP.EQ.6) THEN
          IF (NEWOPT) TOLD = ROPT(IOPT(LP+1))
          GO TO 480

      END IF

      IF (JP.EQ.7) THEN
C     CHANGE TOLERANCE FOR RELATIVE AGREEMENT BETWEEN
C     BEST FUNCTION NORM, LAST FUNCTION NORM AND THE
C     CURRENT FUNCTION NORM.
          IF (NEWOPT) TOLSNR = ROPT(IOPT(LP+1))
          GO TO 480

      END IF

      IF (JP.EQ.8) THEN
C     CHANGE TOLERANCE FOR AGREEMENT BETWEEN PREDICTED
C     VALUE OF RESIDUAL NORM AND THE PREVIOUS VALUE OF
C     THE RESIDUAL NORM.
          IF (NEWOPT) TOLP = ROPT(IOPT(LP+1))
          GO TO 480

      END IF
C     CHANGE THE PRINT LEVEL IN THE ERROR PROCESSOR.
      IF (JP.EQ.9) THEN
          IF (NEWOPT) LEVEL = IOPT(LP+1)
          GO TO 480

      END IF
C     PASS AN OPTION ARRAY TO THE CONSTRAINED LINEAR SOLVER.
C     THIS OPTION IS A POINTER TO THE START OF THE OPTION
C     ARRAY FOR THE SUBPROGRAM.
      IF (JP.EQ.10) THEN
          IF (NEWOPT) IPLS = IOPT(LP+1)
          GO TO 480

      END IF
C     MOVE THE PROCESSING POINTER BY THE VALUE IN THE
C     NEXT ENTRY OF THE OPTION ARRAY.  THIS DEVICE IS
C     INCLUDED SO THAT PASSING OPTIONS TO LOWER LEVEL
C     SUBROUTINES IS EASY TO DO.
      IF (JP.EQ.11) THEN
          IF (NEWOPT) LPDIFF = IOPT(LP+1)
          GO TO 480

      END IF
C     SAW AN OPTION (OR GARBAGE) THAT IS NOT ON THE LIST.
      XMESS =
     .'DQEDIP. INVALID OPTION PROCESSED. I1=IOPT(*) ENTRY. I2=IOPT(I1).'
      NERR = 08
      IGO = 16
      CALL CHRCNT(XMESS,NMESS)
      CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,LP,IOPT(LP),0,RDUM,RDUM)
      IFLAG = 0
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 440
C     END PROCEDURE
  490 CONTINUE
      GO TO IGOPOA

 9001 FORMAT (A4,1P10D12.4/ (4X,10D12.4))
 9011 FORMAT ('0+ITER.=',I3,' FC=',1PD10.4,' PV=',1PD10.4,' K=',I4,
     .  ' KL=',I4,' FB=',1PD10.4/12X,'AL=',1PD14.4,' BB=',1PD14.4)
 9021 FORMAT (' LINEAR RESIDUAL.GE.CURRENT F. QUITTING.',1P2D12.5)
      END
      SUBROUTINE DQEDMN(DQEDEV,MEQUA,NVARS,MCON,IND,BL,BU,X,FJAC,LDFJAC,
     .                  FB,IGO,IOPT,ROPT,IWA,WA,DX,XB,B,BB,BLB,BUB,INDB,
     .                  NPMAX,ZP,XP,QC,MDQC,PJ,WJ,LDWJ,GR,DXL)
C***BEGIN PROLOGUE  DQEDMN
C***REFER TO  DQED
C***ROUTINES CALLED  DGECO,DGESL,IDAMAX,DNRM2,I1MACH,
C                    DQEDEV,DQEDGN,D1MACH,DAXPY,DSCAL,DCOPY,
C                    DVOUT,DDOT,CHRCNT,XERRWV
C***END PROLOGUE  DQEDMN
C     REVISED 861216-1100
C     REVISED YYMMDD-HHMM
      DOUBLE PRECISION BL(*),BU(*),X(*),FJAC(LDFJAC,*)
      DOUBLE PRECISION DX(*),XB(*),B(*),BB(*)
      DOUBLE PRECISION BLB(*),BUB(*),ROPT(*),WA(*),PJ(*)
      DOUBLE PRECISION WJ(LDWJ,*),GR(*),DXL(*)
C     ARRAYS TO HOLD VALUES FOR 2-ND ORDER TERMS .
      DOUBLE PRECISION ZP(MEQUA,NPMAX),XP(NVARS,NPMAX),QC(MDQC,NPMAX)
C     SCALARS:
      DOUBLE PRECISION AJN,ALB,ALFAC,ALPHA,AUB,BBOOST,BOLD,CHG
      DOUBLE PRECISION CHGFAC,COLNRM,COND,COSL,COSM,COSQ,DFN,DXNRM
      DOUBLE PRECISION C1516,FB,FC,FL,PV,PVL,RC,RG,RB,PB,PD
      DOUBLE PRECISION RCOND,SEMIBG,T,TOLD,TOLF,TOLP,TOLSNR,TOLUSE
      DOUBLE PRECISION TOLX,TT,TWO,T2,ZERO,ONE,ZN,SS,SC
      DOUBLE PRECISION D1MACH,DDOT,DNRM2
      INTEGER IND(*),INDB(*),IOPT(*),IWA(*)
      LOGICAL RETREA,TERM,FULNWT,USEQ,NEWBST,PASSB,NOQUAD,NEWOPT
      LOGICAL LINMOD,REVERS,USEQL,MUSTCN,JACTRI
      CHARACTER XMESS*128
      DATA IFLAG/0/
C--PROCEDURES--
C -NAME------TYPE--------ARGS------CLASS-----
C
C  DGECO                    6      SUBROUTINE
C  DGESL                    6      SUBROUTINE
C  IDAMAX    INTEGER        3      FUNCTION
C  DNRM2     REAL           3      FUNCTION
C  DQEDEV                   6      DUMMY-SUBR
C  DQEDGN                  15      SUBROUTINE
C  I1MACH    INTEGER        1      FUNCTION
C  D1MACH    REAL           1      FUNCTION
C  DAXPY                    6      SUBROUTINE
C  DSCAL                    4      SUBROUTINE
C  DCOPY                    5      SUBROUTINE
C  SVOUT                    4      SUBROUTINE
C  DDOT      REAL           5      FUNCTION
C  CHRCNT                   2      SUBROUTINE
C  XERRWV                  10      SUBROUTINE
C DQEDMN:
C GLOSSARY OF VARIABLES. NOTATION:
C DUMMY-ARG A DUMMY ARGUMENT, THAT IS AN ARGUMENT TO THIS PROG. UNIT.
C /S$A$V$E/ SAV DENOTES THAT THIS VARIABLE IS LOCAL TO THE ROUTINE
C               AND IS SAVED BETWEEN CALLS TO IT.
C INTEGER, REAL, DOUBLE PRECISION, LOGICAL, CHARACTER
C               THE TYPES OF THE VARIABLES.
C ADJ-ARR AN ADJUSTABLE ARRAY, THAT IS AN ARGUMENT TO THIS PROG. UNIT.
C NAME      MEMORY STATUS  TYPE     ARGUMENT   USES AND COMMENTS.
C                                    STATUS
C ----      -------------  ----     --------   ------------------
C AJN        /S$A$V$E/ SAV REAL              NORM GRAD VECTOR
C ALB        /S$A$V$E/ SAV REAL              TEMP BOUND VALUE
C ALFAC      /S$A$V$E/ SAV REAL              TRUST REGION FACTOR
C ALPHA      /S$A$V$E/ SAV REAL              TRUST REGION FACTOR
C AUB        /S$A$V$E/ SAV REAL              TEMP BOUND VALUE
C B          DUMMY-ARG     REAL      ADJ-ARY CURRENT TRUST BOUNDS
C BB         DUMMY-ARG     REAL      ADJ-ARY TRUST BOUNDS AT BEST
C BBOOST     /S$A$V$E/ SAV REAL              FACTOR TO BOOST BOUNDS
C BL         DUMMY-ARG     REAL      ADJ-ARY USER LOWER BOUNDS
C BLB        DUMMY-ARG     REAL      ADJ-ARY MODEL LOWER BOUNDS
C BOLD       /S$A$V$E/ SAV REAL              TEMP LAST BOUND VALUE
C BU         DUMMY-ARG     REAL      ADJ-ARY USER UPPER BOUNDS
C BUB        DUMMY-ARG     REAL      ADJ-ARY MODEL UPPER BOUNDS
C CHG        /S$A$V$E/ SAV REAL              TRUST REGION FACTOR
C CHGFAC     /S$A$V$E/ SAV REAL              TRUST REGION FACTOR
C COLNRM     /S$A$V$E/ SAV REAL              TEMP JACOBIAN COL NORM
C COND       /S$A$V$E/ SAV REAL              MAX COND NUMBER, QUAD TERM
C COSL       /S$A$V$E/ SAV REAL              COSINE, GRAD AND LIN STEP
C COSM       /S$A$V$E/ SAV REAL              COSINE, LIN AND QUAD STEP
C COSQ       /S$A$V$E/ SAV REAL              COSINE, GRAD AND QUAD STEP
C C1516      /S$A$V$E/ SAV REAL              THE CONSTANT 15/16
C DFN        /S$A$V$E/ SAV REAL              NORM J*DX/NORM DEL F
C DX         DUMMY-ARG     REAL      ADJ-ARY THE CHANGE, X=X-DX
C DXL        DUMMY-ARG     REAL      ADJ-ARY A CANDIDATE DX, LIN MODEL
C DXNRM      /S$A$V$E/ SAV REAL              NORM OF DX
C FB         DUMMY-ARG     REAL              NORM AT THE BEST X
C FC         /S$A$V$E/ SAV REAL              NORM AT THE CURRENT X
C FJAC       DUMMY-ARG     REAL      ADJ-ARY JACOBIAN ARRAY
C FL         /S$A$V$E/ SAV REAL              NORM AT THE LAST X
C FULNWT     /S$A$V$E/ SAV LOGICAL           FLAG, TOOK FULL STEP
C GR         DUMMY-ARG     REAL      ADJ-ARY GRADIENT VECTOR
C I          /S$A$V$E/ SAV INTEGER           DUMMY LOOP VARIABLE
C ICASE      /S$A$V$E/ SAV INTEGER           GATE VARIABLE
C IFLAG      /S$A$V$E/ SAV INTEGER           INTERNAL FIRST TIME FLAG
C IGO        DUMMY-ARG     INTEGER           DIRECT USER, PROGRAM ACTION
C NERR       /S$A$V$E/ SAV INTEGER           ERROR PROCESSOR NUMBER
C NEWBST     /S$A$V$E/ SAV LOGICAL           FLAG, GOT A NEW BEST
C NEWOPT     /S$A$V$E/ SAV LOGICAL           FLAG, SET AN OPTION
C NIT        /S$A$V$E/ SAV INTEGER           ITERS, QUAD MODEL SOLVING
C NMESS      /S$A$V$E/ SAV INTEGER           LENGTH OF ERROR MESSAGE
C NOQUAD     /S$A$V$E/ SAV LOGICAL           FLAG, SUPPRESS QUAD MODEL
C NOUT       /S$A$V$E/ SAV INTEGER           UNIT NUMBER, SAMPLE OUTPUT
C NP         /S$A$V$E/ SAV INTEGER           POTENTIAL QUAD TERMS + 1
C NPMAX      DUMMY-ARG     INTEGER           MAX NUMBER QUAD TERMS
C NT         /S$A$V$E/ SAV INTEGER           MIN(NVARS+1,MEQUA-1)
C NTTERM     /S$A$V$E/ SAV INTEGER           NUMBER QUAD TERMS USED
C NV         /S$A$V$E/ SAV INTEGER           NUMBER OF VARIABLES, MODEL
C NVARS      DUMMY-ARG     INTEGER           NUMBER OF USER VARIABLES
C ONE        /S$A$V$E/ SAV REAL              THE NUMBER 1.
C PASSB      /S$A$V$E/ SAV LOGICAL           FLAG, USER GAVE TRUST BNDS
C PB         /S$A$V$E/ SAV REAL              PREDICTED RESIDUAL, AT BEST
C PD         /S$A$V$E/ SAV REAL              1.5(FB+PB(PB/FB)) - 4PB
C PJ         DUMMY-ARG     REAL      ADJ-ARY J**T*F AND J*DX AT X
C PV         /S$A$V$E/ SAV REAL              PREDICTED RESIDUAL, CURRENT
C PVL        /S$A$V$E/ SAV REAL              PREDICTED RESIDUAL, LIN MOD
C QC         DUMMY-ARG     REAL      ADJ-ARY QUAD MODEL COEFFICIENTS
C RB         /S$A$V$E/ SAV REAL              ALPHA AFTER NOT A BEST
C RC         /S$A$V$E/ SAV REAL              RETREAT COEFFICIENT
C RCOND      /S$A$V$E/ SAV REAL              RECIPROCAL, QUAD TERM COEFF
C RDUM       /S$A$V$E/ SAV REAL              DUMMY VARIABLE, ERROR PROC
C RETREA     /S$A$V$E/ SAV LOGICAL           FLAG, WILL RETREAT
C REVERS     /S$A$V$E/ SAV LOGICAL           FLAG, REVERSE COMMUNICATION
C IGOELM     /S$A$V$E/ SAV INTEGER           BACK FROM LIN MOD EVALUATE
C IGOEQM     /S$A$V$E/ SAV INTEGER           BACK FROM QUAD MOD EVALUATE
C IGOTFC     /S$A$V$E/ SAV INTEGER           BACK FROM TEST FOR CONVERGE
C IGOW       /S$A$V$E/ SAV INTEGER           DIRECT, SOLVE QUAD MODEL
C IND        DUMMY-ARG     INTEGER   ADJ-ARY USER BOUND INDICATORS
C INDB       DUMMY-ARG     INTEGER   ADJ-ARY INTERNAL BOUND INDICATORS
C IOPT       DUMMY-ARG     INTEGER   ADJ-ARY USER OPTION ARRAY
C IPLS       /S$A$V$E/ SAV INTEGER
C IPRINT     /S$A$V$E/ SAV INTEGER           WANT OUTPUT IF .GT. 0
C ITERS      /S$A$V$E/ SAV INTEGER           NUMBER OF ITERATIONS
C ITMAX      /S$A$V$E/ SAV INTEGER           MAX NUMBER OF ITERATIONS
C IWA        DUMMY-ARG     INTEGER   ADJ-ARY WORKING ARRAY
C J          /S$A$V$E/ SAV INTEGER           DUMMY LOOP VARIABLE
C JACTRI     /S$A$V$E/ SAV LOGICAL           FLAG, JACOBIAN IS TRIANGLE
C JK         /S$A$V$E/ SAV INTEGER           TEMP LOOP VARIABLE
C JP         /S$A$V$E/ SAV INTEGER           OPTION ARRAY POINTER
C K          /S$A$V$E/ SAV INTEGER           DUMMY LOOP VARIABLE
C KL         /S$A$V$E/ SAV INTEGER           DIRECT STATE OF DQED
C KP         /S$A$V$E/ SAV INTEGER           OPTION ARRAY POINTER
C L          /S$A$V$E/ SAV INTEGER           DUMMY LOOP VARIABLE
C LDFJAC     DUMMY-ARG     INTEGER           ROW DIMENSION OF FJAC
C LDWJ       DUMMY-ARG     INTEGER           ROW DIMENSION OF WJ
C LEVEL      /S$A$V$E/ SAV INTEGER           ERROR PROC RESPONSE LEVEL
C LINMOD     /S$A$V$E/ SAV LOGICAL           FLAG, SOLVING LINEAR PROBLEM
C LK         /S$A$V$E/ SAV INTEGER           MIN(MEQUA,NVARS+1)
C LP         /S$A$V$E/ SAV INTEGER           OPTION ARRAY POINTER
C LPDIFF     /S$A$V$E/ SAV INTEGER           OPTION ARRAY INCREMENT
C MCON       DUMMY-ARG     INTEGER           NUMBER, LINEAR CONSTRAINTS
C MCONST     /S$A$V$E/ SAV INTEGER           NUMBER, MODEL CONSTRAINTS
C MDQC       DUMMY-ARG     INTEGER           ROW DIM OF QC(,)
C ME         /S$A$V$E/ SAV INTEGER           NUMBER, MODEL EQUATIONS
C MEQUA      DUMMY-ARG     INTEGER           NUMBER, USER EQUATIONS
C MK         /S$A$V$E/ SAV INTEGER           0 OR MIN(MEQUA,NVARS+NP+1)
C MUSTCN     /S$A$V$E/ SAV LOGICAL           FLAG, MUST TAKE FULL STEP
C NALL       /S$A$V$E/ SAV INTEGER           MCON+NVARS
C RG         /S$A$V$E/ SAV REAL              RETREAT GAUGE
C ROPT       DUMMY-ARG     REAL      ADJ-ARY USER PASSED OPTION DATA
C SEMIBG     /S$A$V$E/ SAV REAL              CONSTANT 1.D+10
C T          /S$A$V$E/ SAV REAL              TEMP VARIABLE
C TERM       /S$A$V$E/ SAV LOGICAL           FLAG, STOP THIS PROBLEM
C TOLD       /S$A$V$E/ SAV REAL              REL TOLERANCE ON CHANGE
C TOLF       /S$A$V$E/ SAV REAL              TOLERANCE ON FUNCTION
C TOLP       /S$A$V$E/ SAV REAL              TOLERANCE FOR MIN FLAG
C TOLSNR     /S$A$V$E/ SAV REAL              TOLERANCE FOR MIN FLAG
C TOLUSE     /S$A$V$E/ SAV REAL              TOLERANCE FOR QUAD MODEL
C TOLX       /S$A$V$E/ SAV REAL              ABS TOLERANCE ON CHANGE
C TT         /S$A$V$E/ SAV REAL              TEMP VARIABLE
C TWO        /S$A$V$E/ SAV REAL              CONSTANT 2.
C T2         /S$A$V$E/ SAV REAL              REL STEP WITHIN TRUST REG
C USEQ       /S$A$V$E/ SAV LOGICAL           USE QUAD MODEL THIS STEP
C USEQL      /S$A$V$E/ SAV LOGICAL           USED QUAD MODEL LAST STEP
C UT         DUMMY-ARG     REAL      ADJ-ARY WORKING ARRAY
C WA         DUMMY-ARG     REAL      ADJ-ARY WORKING ARRAY
C WJ         DUMMY-ARG     REAL      ADJ-ARY MODEL JACOBIAN ARRAY
C X          DUMMY-ARG     REAL      ADJ-ARY SOLUTION ARRAY
C XB         DUMMY-ARG     REAL      ADJ-ARY BEST VALUES OF X
C XMESS      /S$A$V$E/ SAV CHAR*128          TEMP FOR ERROR MESSAGE
C XP         DUMMY-ARG     REAL      ADJ-ARY WORKING ARRAY
C ZERO       /S$A$V$E/ SAV REAL              CONSTANT 0.
C ZN         /S$A$V$E/ SAV REAL              NORM J*DX/NORM DF
C ZP         DUMMY-ARG     REAL      ADJ-ARY WORKING ARRAY
C
      IF (IFLAG.NE.0) GO TO 50

      LK = MIN(MEQUA,NVARS+1)
      NT = MIN(NVARS+1,MEQUA-1)
      ZERO = 0.D0
      ONE = 1.D0
      TWO = 2.D0
C     DO(PROCESS OPTION ARRAY)
      GO TO 1100

   10 CONTINUE
C     DO(INITIALIZE OTHER VALUES)
      GO TO 1030

   20 CONTINUE
C     SET SO X(*)-DX(*) UPDATE WILL BE CORRECT FIRST TIME.
      DX(1) = ZERO
      CALL DCOPY(NVARS,DX,0,DX,1)
      K = 0
C     D1MACH(2)="INFINITY" ON THIS MACHINE.
      FB = D1MACH(2)
      DXNRM = FB
      FL = ZERO
C
C     MODEL PROBLEM RESIDUAL.
      PV = ZERO
      PVL = ZERO
      RETREA = .FALSE.
      FULNWT = .FALSE.
      TERM = .FALSE.
C     DO FOREVER
   30 CONTINUE
      ITERS = ITERS + 1
      IF (RETREA) THEN
C     MUST RETREAT TO BEST X VALUES.
          CALL DCOPY(NVARS,XB,1,X,1)
          K = 0
          KL = -1
          FL = FB

      ELSE
          KL = K
          DO 40 J = 1,NVARS
             X(J) = X(J) - DX(J)
   40     CONTINUE
      END IF

      IF (TERM) THEN
          IFLAG = 0
C     EXIT FOREVER
          GO TO 840

      END IF

      IFLAG = 1
      IGO = 1
      IF (NP.EQ.NPMAX-1 .AND. NP.LT.NVARS) IGO = -1
      IF (REVERS) THEN
C     THERE ARE TWO POSSIBLE WAYS TO GET FUNCTION AND DERIVATIVE
C     VALUES FROM THE USER.  THE OPTIONAL WAY IS REVERSE COMMUNICATION.
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 1020

      ELSE
C     THE NOMINAL WAY IS THROUGH FORWARD COMMUNICATION.
          CALL DQEDEV(X,FJAC,LDFJAC,IGO,IOPT,ROPT)
      END IF

   50 CONTINUE
C     IF IGO HAS BEEN CHANGED BY THE USER TO A VALUE .GT. 1, THEN
C     THIS IS AN ABORT SIGNAL.  STOP UNLESS IT = 99.
      IF (IGO.EQ.99) THEN
C     IF IGO = 99 THE EVALUATION CAN'T BE PERFORMED.
C     WE FORCE A RETREAT AND RESTART IN THIS CASE.
          DO 70 I = MCON + 1, MCON + MEQUA
             FJAC(I,NVARS+1) = FC
             DO 60 J = 1,NVARS
                FJAC(I,J) = ZERO
   60        CONTINUE
   70     CONTINUE
C     A RETREAT IS FORCED TO OCCUR WITH THIS ASSIGNMENT.
          RETREA = .TRUE.

          END IF

          FC = DNRM2(MEQUA,FJAC(MCON+1,NVARS+1),1)

      IF (IGO.GT.1 .AND. IGO.NE.99) THEN
          IFLAG = 0
          CALL DCOPY(NVARS,XB,1,X,1)
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 1020

      END IF
C     SAVE PAST FUNCTION AND VARIABLE VALUES.
C     DO NOT UPDATE THE PAST POINTS UNLESS THERE IS A
C     SIGNIFICANT CHANGE IN THE X(*) VALUES.
      IF (NP.GE.0) THEN
          IF (DXNRM.GT.TOLUSE*DNRM2(NVARS,X,1)) THEN
              LP = NVARS
              IF ( .NOT. NOQUAD) NP = MIN(NP,NPMAX-1,LP) + 1
              DO 150 J = NP - 1,1,-1
C     SAVE THE PAST VALUES OF THE VARIABLES.
C     DIFFERENCES ARE LATER COMPUTED USING THESE VALUES.
                 CALL DCOPY(NVARS,XP(1,J),1,XP(1,J+1),1)
C     SAVE THE PAST FUNCTION VALUES IN ZP( , ).
C     DIFFERENCES ARE LATER COMPUTED USING THESE VALUES.
                 CALL DCOPY(MEQUA,ZP(1,J),1,ZP(1,J+1),1)
  150         CONTINUE
          END IF

      END IF
C     PUT IN THE PRESENT VALUES OF THE VARIABLES.
      CALL DCOPY(NVARS,X,1,XP(1,1),1)
C     PUT IN THE PRESENT VALUES OF THE FUNCTIONS.
      CALL DCOPY(MEQUA,FJAC(MCON+1,NVARS+1),1,ZP(1,1),1)
C     THIS STATEMENT HAS THE EFFECT OF A FIRST TIME FLAG.
      NP = MAX(NP,0)
C     COMPUTE THE COSINES OF THE PAST MOVES WITH THE MOST CURRENT MOVE.
      DO 170 L = 2,NP
         DO 160 J = 1,NVARS
            QC(J,L) = XP(J,L) - XP(J,1)
  160    CONTINUE
  170 CONTINUE
      L = 3
  180 CONTINUE
C     DO WHILE(L.LE.NP)
      IF ( .NOT. (L.LE.NP)) GO TO 200
C     CALCULATE THE DIRECTION COSINES OF THE PAST MOVES.
      T = DDOT(NVARS,QC(1,2),1,QC(1,L),1)
      TT = DNRM2(NVARS,QC(1,2),1)*DNRM2(NVARS,QC(1,L),1)
      IF (TT.GT.ZERO) THEN
          T = T/TT

      ELSE
          T = ONE
      END IF

      IF (IPRINT.GT.0) THEN
          WRITE (NOUT,
     .      '('' PAST MOVE NUMBER, COSINE OF MOVE'',I3,2X,F6.2)') L - 2,
     .      T
      END IF

      IF(ABS(T) .GT. .98) THEN
C     DISCARD PAST INFORMATION ASSOCIATED WITH THIS MOVE IF CLOSE TO
C     A PAST MOVE.
          DO 190 J = L,NP - 1
             CALL DCOPY(MEQUA,ZP(1,J+1),1,ZP(1,J),1)
             CALL DCOPY(NVARS,XP(1,J+1),1,XP(1,J),1)
             CALL DCOPY(NVARS,QC(1,J+1),1,QC(1,J),1)
  190     CONTINUE
          NP = NP - 1
C     CYCLE WHILE
          GO TO 180

      END IF

      L = L + 1
C     END WHILE
      GO TO 180

  200 CONTINUE
C     COMPUTE FUNCTION DIFFERENCES IN QC( , ).
      DO 220 J = 1,NP - 1
         DO 210 I = 1,MEQUA
            QC(I,J+1) = ZP(I,J+1) - ZP(I,1)
  210    CONTINUE
  220 CONTINUE
C     NOW HAVE F(PAST)-F(CURRENT) IN QC( , ), COLS. 2,...,NP USED.
C     COMPUTE NORM OF DIFFERENCE OF FUNCTION VALUES.
      IF (NP.GT.1) THEN
          DFN = DNRM2(MEQUA,QC(1,2),1)

      ELSE
          DFN = ZERO
      END IF

      DO 240 I = 1,NP - 1
         DO 230 J = 1,NVARS
C     NEXT ADD PRODUCT OF JACOBIAN AND PAST VARIABLE DIFFERENCES.
            CALL DAXPY(MEQUA,- (XP(J,I+1)-XP(J,1)),FJAC(MCON+1,J),1,
     .                 QC(1,I+1),1)
  230    CONTINUE
  240 CONTINUE
C     DO FOREVER
  250 CONTINUE
C     COMPUTE THE SYMMETRIC MATRIX WHOSE ENTRIES ARE THE
C     SQUARES OF THE DOT PRODUCTS OF THE PAST DIRECTIONS.
C     THIS MATRIX IS REQUIRED TO OBTAIN THE QUADRATIC TERMS
C     ASSOCIATED WITH INTERPOLATING TO PAST FUNCTION VALUES.
      DO 280 L = 2,NP
         DO 270 J = L,NP
            T = ZERO
            DO 260 I = 1,NVARS
               T = T + (XP(I,J)-XP(I,1))* (XP(I,L)-XP(I,1))
  260       CONTINUE
            WJ(L-1,J-1) = T
            WJ(J-1,L-1) = T
  270    CONTINUE
  280 CONTINUE
C     COMPUTE NORM OF REMAINDER INCLUDING LINEAR TERMS,
C     USING THE LAST MOVE.
      USEQ = NP .GT. 1 .AND. .NOT. RETREA
      ZN = ONE
      IF (NP.GT.1) THEN
          ZN = DNRM2(MEQUA,QC(1,2),1)
C     COMPUTE RATIO OF Z TERMS TO CURRENT F VALUE..
          IF (USEQ) USEQ = (ZN .GT.1.D-4*DFN .AND. ZN .LT. DFN*.75D0)
     .              .OR. USEQL
          IF (DFN.GT.ZERO) THEN
              ZN = ZN/DFN

          ELSE
              ZN = ONE
          END IF

          IF (IPRINT.GT.0) THEN
              CALL DVOUT(1,ZN,'('' RATIO OF Z TERM TO PAST DF NORM'')',
     .                   4)
          END IF
C     SCALE THE MATRIX (MATRIX := D*MATRIX*D, WHERE D**2 = RECIPROCAL
C     OF THE DIAGONAL TERMS OF THE MATRIX.
          DO 290 I = 1,NP - 1
             DXL(I) = WJ(I,I)
             IF (DXL(I).EQ.ZERO) THEN
                 NP = I
                 GO TO 250

             ELSE
                 DXL(I) = ONE/DXL(I)
             END IF

  290     CONTINUE
          DO 310 I = 1,NP - 1
             DO 300 J = 1,NP - 1
                WJ(I,J) = (WJ(I,J)*DXL(I))* (WJ(I,J)*DXL(J))
  300        CONTINUE
  310     CONTINUE
C     USE THE LINPACK ROUTINES DGECO(), DGESL() TO OBTAIN
C     THE COEFFICIENTS OF THE QUADRATIC TERMS, ONE ROW AT
C     A TIME.
          CALL DGECO(WJ,LDWJ,NP-1,IWA,RCOND,WA)
          IF (IPRINT.GT.0) WRITE (NOUT,
     .        '('' RCOND FROM DGECO() = '',2X,       1PD15.4)') RCOND
          IF (COND*RCOND.LT.ONE) THEN
              NP = NP - 1
C     CYCLE FOREVER
              GO TO 250

          END IF

          DO 340 I = 1,MEQUA
C     COPY A ROW OF THE INTERPOLATED DATA TO A WORKING ARRAY.
C     USE THIS ARRAY TO OBTAIN A ROW OF THE QUADRATIC TERMS.
             CALL DCOPY(NP-1,QC(I,2),MDQC,WA,1)
C     SCALE THE RIGHT HAND SIDE DATA.
             DO 320 J = 1,NP - 1
                WA(J) = WA(J)*DXL(J)
  320        CONTINUE
             CALL DGESL(WJ,LDWJ,NP-1,IWA,WA,0)
C     RESCALE THE SOLUTION DATA.
             DO 330 J = 1,NP - 1
                WA(J) = WA(J)*DXL(J)
  330        CONTINUE
C     THE FOLLOWING SIGN CHANGE COMES FROM A CHANGE
C     OF SIGN IN THE INNER LOOP MODEL PROBLEM.
             CALL DSCAL(NP-1,-TWO,WA,1)
             CALL DCOPY(NP-1,WA,1,QC(I,2),MDQC)
  340     CONTINUE
      END IF

  350 CONTINUE
C     EXIT FOREVER
C     END FOREVER
C     NOW HAVE THE QUADRATIC TERMS COMPUTED.
C     NEXT WILL TRIANGULARIZE THE JACOBIAN TO SAVE SPACE
C     WHEN USING THE QUADRATIC MODEL.
      IF (JACTRI) THEN
C     CONSTRUCT AND THEN APPLY PLANE ROTATIONS
C     TO ACHIEVE UPPER TRIANGULAR FORM.  THIS LOOP
C     AFFECTS THE JACOBIAN AND RIGHT HAND SIDE.
          DO 370 J = 1,NT
             DO 360 I = J + 1,MEQUA
                CALL DROTG(FJAC(MCON+J,J),FJAC(MCON+I,J),SC,SS)
                CALL DROT(NVARS-J+1,FJAC(MCON+J,J+1),LDFJAC,
     .                              FJAC(MCON+I,J+1),LDFJAC,SC,SS)
C     NOW APPLY THE TRANSFORMATION TO THE QUADRATIC TERMS.
                CALL DROT(NP-1,QC(J,2),MDQC,
     .                         QC(I,2),MDQC,SC,SS)

                FJAC(MCON+I,J) = ZERO
  360        CONTINUE
  370     CONTINUE
C     NOW WE FINISH TRIANGULARIZING THE QUADRATIC TERMS.
C     NOTE THAT THIS DOES NOT AFFECT THE RIGHT HAND SIDE.
          DO 390 L = 1,NP - 1
             DO 395 I=NVARS+L+2,MEQUA
             CALL DROTG(QC(NVARS+L+1,L+1),QC(I,L+1),SC,SS)
             CALL DROT(NP-L-1,QC(NVARS+L+1,MIN(L+2,NPMAX)),MDQC,
     .                        QC(I,MIN(L+2,NPMAX)),MDQC,SC,SS)
             QC(I,L+1) = ZERO
  395 CONTINUE
  390     CONTINUE
      END IF
C     COMPUTE CURRENT NORM OF J**T*F(X).
      DO 400 J = 1,NVARS
         IF (JACTRI) THEN
             JK = J

         ELSE
             JK = MEQUA
         END IF

         PJ(J) = DDOT(JK,FJAC(MCON+1,J),1,FJAC(MCON+1,NVARS+1),1)
  400 CONTINUE
      AJN = DNRM2(NVARS,PJ,1)
C SAVE J**T*F FOR DIRECTION TESTING WITH LINEAR AND QUADRATIC MOVES.
      IF (AJN.GT.ZERO) CALL DSCAL(NVARS,ONE/AJN,PJ,1)
      CALL DCOPY(NVARS,PJ,1,GR,1)
      NEWBST = FC .LT. FB
      IF (NEWBST) K = 0
      IF (K.EQ.0) THEN
          PB = ZERO
          PD = D1MACH(2)
C     WANT TO POSITION AT BEST X VALUES.
          IF ( .NOT. RETREA) THEN
              FB = FC
          END IF

          GO TO (410,430,450),2 - KL

          GO TO 470
C     CASE 1
  410     CONTINUE
C     IMMEDIATELY GOT A NEW BEST X.
          RG = ZERO
          IF (T2.LE.0.25D0) THEN
              BBOOST = ONE
              CHG = MAX(4.D0*T2,.1D0)
          END IF

          DO 420 J = 1,NVARS
             BB(J) = CHG*BB(J)
  420     CONTINUE
          T = .25D0/ (ALFAC-1.D0)
          ALPHA = (ZN+ALFAC*T)/ (ZN+T)
          ALFAC = 1.5*ALPHA
          BBOOST = MIN(1.5D0*ALPHA*BBOOST,SEMIBG)
          GO TO 490
C     CASE 2
  430     CONTINUE
          USEQL = .FALSE.
          RG = ZERO
C     AT THE INITIAL X.
          ALFAC = 256.D0
          DO 440 J = 1,NVARS
             IF ( .NOT. PASSB) THEN
                 BB(J) = -X(J)
             END IF

             IF (BB(J).EQ.ZERO) THEN
                 IF (JACTRI) THEN
                     JK = J

                 ELSE
                     JK = MEQUA
                 END IF

                 COLNRM = DNRM2(JK,FJAC(MCON+1,J),1)
                 IF (COLNRM.NE.ZERO) THEN
                     BB(J) = DDOT(JK,FJAC(MCON+1,J),1,
     .                       FJAC(MCON+1,NVARS+1),1)
                     BB(J) = -MAX(ABS(BB(J))/COLNRM/COLNRM,FC/COLNRM)

                 ELSE
                     BB(J) = -ONE
                 END IF

             END IF

             XB(J) = X(J)
             B(J) = BB(J)
  440     CONTINUE
          ALPHA = ONE
          BBOOST = 0.5D0
C     EXIT IF
          GO TO 520
C     CASE 3
  450     CONTINUE
C     RETREAT TO BEST X.
          IF (ALFAC.NE.256.D0) THEN
              ALPHA = MIN(4.D0/ALFAC,.25D0)
              ALFAC = 1.25D0

          ELSE
              ALPHA = .25D0*ALPHA
          END IF

          BBOOST = .25D0
          USEQL = .FALSE.
          DO 460 J = 1,NVARS
C
C     THE NEXT LINES HAVE THE EFFECT OF REDUCING THE BOUNDS
C     AT THE CURRENT BEST TO ABOUT 1/16 THEIR CURRENT VALUES
C     PROVIDED THE CURRENT BOUNDS ARE RELATIVELY SMALL.  IF
C     THE CURRENT BOUNDS ARE RELATIVELY LARGE, THE BOUNDS AT
C     THE BEST ARE LEFT ABOUT THE SAME.
             T = ABS(B(J))
             TT = ABS(BB(J))
             T = (T+.25D0*TT)/ (T+4.D0*TT)
             B(J) = T*BB(J)
             BB(J) = B(J)
             DX(J) = ZERO
  460     CONTINUE
C     EXIT IF
          GO TO 520
C     CASE OTHER
  470     CONTINUE
C     NOT IMMEDIATELY A BEST X.
          RB = ZERO
          DO 480 J = 1,NVARS
             RB = .125D0*MAX(RB,ABS((XB(J)-X(J))/BB(J)))
  480     CONTINUE
          ALPHA = RB
          ALFAC = TWO
          BBOOST = (1.D0+RG)/ (.25D0+2.D0*RG)
C     END CASE
  490     CONTINUE
          DO 500 J = 1,NVARS
             DX(J) = XB(J) - X(J)
             IF (DX(J).EQ.ZERO) THEN
                 B(J) = ALPHA*BB(J)

             ELSE
                 XB(J) = X(J)
                 B(J) = SIGN(ALPHA*BB(J),DX(J)) + BBOOST*DX(J)
             END IF

             BB(J) = SIGN(MIN(SQRT(D1MACH(2)),ABS(B(J))),B(J))
  500     CONTINUE

      ELSE
C     MUST MAKE SURE THAT PD GETS SET TO A REASONABLE VALUE.
C     COMPUTE A GAUGE FOR RETREATING IF DID NOT
C     GET A NEW BEST.
          ALPHA = (.5D0*ZN+1.D0)/ (ZN+1.D0)
          CHG = ALPHA*CHG
          IF (K.EQ.1) THEN
              CHG = MIN(CHG,T2)
              CHG = MAX(CHG,.1D0)
              PB = PV
              PD = 1.5D0* (FB+PB* (PB/FB)) - 4.D0*PB
          END IF

          DO 510 J = 1,NVARS
             B(J) = CHG*B(J)
             IF (K.EQ.1) BB(J) = B(J)
  510     CONTINUE
      END IF

  520 CONTINUE
C     DO(TEST FOR CONVERGENCE)
      ASSIGN 530 TO IGOTFC
      GO TO 980

  530 CONTINUE
      IF (TERM) THEN
          IFLAG = 0
C     EXIT FOREVER
          GO TO 840

      END IF

      K = K + 1
C     SOLVE MODEL BOUNDED PROBLEM.
      DO 590 J = 1,NVARS
         IF (B(J).LT.ZERO) THEN
             ALB = B(J)
             IF (DX(J).EQ.ZERO) THEN
C     THIS CASE IS REQD. TO AVOID USING BUB(*) AT THE INITIAL PT.
                 AUB = -C1516*ALB

             ELSE
                 AUB = MIN(-C1516*ALB,-DX(J)+BUB(J))
             END IF

         ELSE
             AUB = B(J)
             IF (DX(J).EQ.ZERO) THEN
                 ALB = -C1516*AUB

             ELSE
                 ALB = MAX(-C1516*AUB,-DX(J)+BLB(J))
             END IF

         END IF
c
c     THIS NEXT CODE, ENDING WITH ***, POINTS THE BOX TOWARDS THE BEST
c     VALUE OF X WHEN NOT AT A NEW BEST.
c
         IF (K.GE.2) THEN
             IF (XB(J).GT.X(J)) THEN
                 BUB(J) = BUB(J)*.25D0
                 IF (X(J)-BLB(J).GT.XB(J)) THEN
                     BLB(J) = BLB(J)*.75D0

                 ELSE
                     BLB(J) = MIN(BLB(J),.75D0* (X(J)-XB(J)))
                 END IF

             ELSE
                 BLB(J) = BLB(J)*.25D0
                 IF (X(J)-BUB(J).LT.XB(J)) THEN
                     BUB(J) = BUB(J)*.75D0

                 ELSE
                     BUB(J) = MAX(BUB(J),.75D0* (X(J)-XB(J)))
                 END IF

             END IF

         END IF
C***
C     RESTRICT THE STEP FURTHER IF USER GIVES BOUNDS.
         ICASE = IND(J)
         GO TO (540,550,560,570),ICASE
C     CASE 1
  540    AUB = MIN(AUB,X(J)-BL(J))
         GO TO 580
C     CASE 2
  550    ALB = MAX(ALB,X(J)-BU(J))
         GO TO 580
C     CASE 3
  560    AUB = MIN(AUB,X(J)-BL(J))
         ALB = MAX(ALB,X(J)-BU(J))
         GO TO 580
C     CASE 4
  570    CONTINUE
C     END CASE
  580    BLB(J) = ALB
C     THIS NEXT LINE IS TO GUARANTEE THAT THE LOWER BOUND
C     IS .LE. THE UPPER BOUND.
         AUB = MAX(AUB,ALB)
         BUB(J) = AUB
         INDB(J) = 3
  590 CONTINUE
C     COMPUTE JACOBIAN*DX AND COMPARE NORM WITH CURRENT FUNCTION.
      IF (NP.GT.1) THEN
          PJ(1) = ZERO
          CALL DCOPY(LK,PJ,0,PJ,1)
          DO 600 J = 1,NVARS
             IF (JACTRI) THEN
                 JK = J

             ELSE
                 JK = MEQUA
             END IF

             CALL DAXPY(JK,DX(J),FJAC(MCON+1,J),1,PJ,1)
  600     CONTINUE
          T = DNRM2(LK,PJ,1)
C   THIS TEST SAYS TO USE THE QUADRATIC MODEL IF
C   THE LAST STEP IS APPROXIMATELY IN THE NULL SPACE
C   OF THE JACOBIAN.
          USEQ = USEQ .OR. T .LT. DFN*0.75D0
          IF (DFN.GT.ZERO) THEN
              DFN = T/DFN

          ELSE
              DFN = ZERO
          END IF

      END IF

      IF (IPRINT.GT.0) THEN
          CALL DVOUT(1,DFN,'('' RATIO OF J*DX NORM TO PAST DF NORM'')',
     .               4)
      END IF
C     CHECK IF QUAD. MODEL IS BEING SUPPRESSED.
      USEQ = USEQ .AND. .NOT. NOQUAD
C     START THE PROCESS USING THE LINEAR MODEL.
      LINMOD = .TRUE.
      MCONST = MCON
CA:   DO FOREVER
  610 CONTINUE
C     COMPUTE THE REQUIRED DIMENSIONS OF THE MODEL PROBLEM.
      IF (LINMOD) THEN
          MK = 0
          ME = MIN(MEQUA,NVARS+1)
C     SET THE INITIAL VALUES FOR THE LINEAR MODEL PROBLEM.
          DX(1) = ZERO
          CALL DCOPY(NVARS,DX,0,DX,1)

      ELSE IF (USEQ) THEN
          MK = MIN(MEQUA,NVARS+NP+1)
          ME = NVARS + MK
          CALL DCOPY(MK,FJAC(MCON+1,NVARS+1),1,DX(NVARS+1),1)

      ELSE
C     EXIT FOREVER(A)
          GO TO 730

      END IF

      NV = NVARS + MK
C     NOTE THAT THE RESIDUALS ARE FREE VARIABLES.
      DO 620 I = NVARS + 1,NV
         INDB(I) = 4
  620 CONTINUE
      NIT = 0
C     THE JACOBIAN, RIGHT SIDE, QUAD. TERMS ARE AN
C     UPPER TRAPAZOIDAL DATA ARRAY.  THIS WILL MAKE SOLVING
C     THE MODEL PROBLEM MORE EFFICIENT.
C     DO FOREVER
  630 CONTINUE
C     CALL A SIMPLIFIED VERSION OF THE ALGORITHM TO SOLVE
C     THE MODEL PROBLEM.  THE FUNCTION AND JACOBIAN
C     ARE COMPUTED FOR THIS SUBPROBLEM DURING THE REVERSE
C     COMMUNICATION REQUESTS.
      CALL DQEDGN(ME,NV,MCONST,INDB,BLB,BUB,DX,WJ,LDWJ,PV,IGOW,
     .            IOPT(IPLS),ROPT,IWA,WA)
C     CHECK FOR AN ERROR THAT WAS SEEN IN THE LOW-LEVEL
C     NONLINEAR SOLVER.
      IF (IGOW.GT.7) THEN
          IGO = IGOW
          IFLAG = 0
C     DO(RETURN TO USER PROGRAM UNIT)
          GO TO 1020

      END IF
C     CLEAR OUT THE WJ(*,*) ARRAY THAT HOLDS
C     THE JACOBIAN FOR THE INNER LOOP PROBLEM.
      WJ(1,1) = ZERO
      CALL DCOPY(LDWJ* (NV+1),WJ,0,WJ,1)
      IF (USEQ .AND. .NOT. LINMOD) WJ(MCONST+1,NVARS+1) = ONE
C     PUT IN A UNIT MATRIX FOR THE PARTIALS
C     WITH RESPECT TO THE RESIDUALS.
      CALL DCOPY(MK,WJ(MCONST+1,NVARS+1),0,WJ(MCONST+1,NVARS+1),LDWJ+1)
      DO 640 I = 1,MCON
C     THE FORM OF THE UPDATE BEING COMPUTED IS X(*)-DX(*).
C     THE VALUE OF DX IS ITSELF COMPUTED AS THE SOLUTION
C     TO A NONLINEAR PROBLEM WITH UPDATES OF THE FORM
C     DX(*)-D(DX)(*).
         CALL DCOPY(NVARS,FJAC(I,1),LDFJAC,WJ(I,1),LDWJ)
         WJ(I,NV+1) = FJAC(I,NVARS+1) - DDOT(NVARS,DX,1,WJ(I,1),LDWJ)
  640 CONTINUE
C     SEE IF USER HAS GIVEN GENERAL CONSTRAINTS.
C     USE THESE CONSTRAINTS TO PLACE EQUIVALENT CONSTRAINTS
C     ON THE CHANGES BEING COMPUTED.
      DO 650 J = 1,MCON
         BLB(NV+J) = BL(J+NVARS)
         BUB(NV+J) = BU(J+NVARS)
         INDB(NV+J) = IND(J+NVARS)
  650 CONTINUE
C     DO(EVALUATE LINEAR MODEL)
      ASSIGN 660 TO IGOELM
      GO TO 940

  660 CONTINUE
C     DO(EVALUATE QUADRATIC MODEL)
      ASSIGN 670 TO IGOEQM
      GO TO 850

  670 CONTINUE
      IF (IGOW.GT.1) THEN
          PV = DNRM2(ME,WJ(MCONST+1,NV+1),1)
          IF (LINMOD) THEN
              PVL = PV
              CALL DCOPY(NVARS,DX,1,DXL,1)
              LINMOD = .FALSE.
C     CYCLE FOREVER(A)
              GO TO 610
C     IF THE PREDICTED NORM IS GREATER THAN THE CURRENT
C     RESIDUAL NORM, DROP THE QUADRATIC MODEL AND USE THE
C     LINEAR MODEL.
          ELSE IF ((PV.GE.FC) .AND. USEQ) THEN
              IF (IPRINT.GT.0) WRITE (NOUT,
     .            '('' ABANDON QUAD. MODEL.'')')
              USEQ = .FALSE.
          END IF

C     EXIT FOREVER(A)
          GO TO 730

      END IF
C     FOR EITHER CASE TRANSFER THE JACOBIAN FOR THE MODEL
C     PROBLEM.  THE TRANSPOSE OF THIS MATRIX IS THE PARTIALS
C     WITH RESPECT TO THE RESIDUALS.
      DO 680 J = 1,NVARS
         CALL DCOPY(MK,WJ(MCONST+1,J),1,WJ(MCONST+MK+J,NVARS+1),LDWJ)
  680 CONTINUE
C     NOW UPDATE THE RESIDUALS FOR BOTH SETS OF MODEL EQUATIONS.
C     IN ROWS 1,...,MK THIS INVOLVES ADDING DX(NVARS+I) TO ROW I.
C     FOR ROWS MK+1,...,ME THIS REQUIRES ADDING MULTIPLES OF THE
C     COLS. OF THE TRANSPOSED JACOBIAN.
      DO 690 I = 1,MK
         T = DX(NVARS+I)
         WJ(MCONST+I,NV+1) = WJ(MCONST+I,NV+1) + T
  690 CONTINUE
C     SYMMETRIZE THE SECOND DERIVATIVE MATRIX.  THIS
C     IS NOT REQUIRED WHEN THE MODEL IS LINEAR, BUT IT
C     DOES NOT HURT THEN.
      IF (USEQ .AND. .NOT. LINMOD) THEN
          DO 710 J = 1,NVARS
             DO 700 I = J,NVARS
                WJ(MCONST+MK+I,J) = WJ(MCONST+MK+J,I)
  700        CONTINUE
  710     CONTINUE
      END IF
C     COMPUTE RESIDUALS ON THE EQUATIONS K*R = 0.
      DO 720 J = NVARS + 1,NV
         CALL DAXPY(NVARS,DX(J),WJ(MCONST+MK+1,J),1,
     .              WJ(MCONST+MK+1,NV+1),1)
  720 CONTINUE
      NIT = NIT + 1
      GO TO 630
C     END FOREVER
C     END FOREVER (A)
  730 CONTINUE
C     COMPUTE THE ANGLES BETWEEN THE LINEAR AND QUADRATIC STEP.
C     TAKE THE ONE, IF THERE IS A CHOICE, CLOSEST TO THE GRADIENT.
C     IF THE QUADRATIC MOVE IS QUITE CLOSE TO THE GRADIENT, TAKE
C     THAT MOVE IN PREFERENCE TO THE LINEAR MOVE.
      COSL = DDOT(NVARS,GR,1,DXL,1)
      T = DNRM2(NVARS,DXL,1)
      IF (T.GT.ZERO) COSL = COSL/T
      COSQ = -ONE
      COSM = -ONE
      TT = ZERO
      IF (USEQ) THEN
          COSQ = DDOT(NVARS,GR,1,DX,1)
          TT = DNRM2(NVARS,DX,1)
          IF (TT.GT.ZERO) COSQ = COSQ/TT
C     COMPUTE THE COSINE OF THE ANGLE BETWEEN THE QUAD. AND
C     LINEAR MOVES.
          IF (T.GT.ZERO .AND. TT.GT.ZERO) COSM = DDOT(NVARS,DX,1,DXL,1)/
     .        T/TT
      END IF

      IF (IPRINT.GT.0) WRITE (NOUT,
     .'('' COS OF QUAD. MOVE AND GRAD., COS OF LIN. MOVE AND GRAD., COS
     . OF EACH MOVE'',/,
     .'' FLAG FOR TRYING QUAD. MOVE.'',/,
     .1P3D12.4,L6)') COSQ,COSL,COSM,USEQ
      IF (IPRINT.GT.0) WRITE (NOUT,
     .'('' LENGTH OF QUAD., THEN LINEAR MOVES'',/,
     .1P2D12.4)') TT,T
C     CHOOSE MOVE PARTIALLY BASED ON ANGLE MOVES MAKE WITH EACH OTHER.
      USEQ = USEQ .AND. (COSM.GT.ZERO .OR. COSL.LT.COSM) .AND. COSQ .GT.
     .       ZERO
      USEQL = USEQ

      IF ( .NOT. USEQ) THEN
          PV = PVL
          CALL DCOPY(NVARS,DXL,1,DX,1)
          NTTERM = 0

      ELSE
          NTTERM = NP - 1
      END IF
C     TEST FOR NOISE IN MODEL PROBLEM SOLN.
      TERM = (PV.GE.FC) .AND. .NOT. RETREA .AND. .NOT. USEQ
      TERM = TERM .AND. MCON .EQ. 0
      TERM = .FALSE.
      IF (TERM) THEN
          IF (IPRINT.GT.0) THEN
              WRITE (NOUT,9021) PV,FC
          END IF
C     VALUE MEANS MODEL RES. .GE. NONLINEAR FUNCTION VALUE.
          IGO = 5
          IFLAG = 0
C     EXIT FOREVER
          GO TO 840

      END IF

      RC = ZERO
      IF (PV.GT.PB) RC = 4.D0* (PV-PB)/ (PV+PD)
C
C     IF USING A QUADRATIC MODEL AND RETREATING SEEMS TO BE
C     NECESSARY, SEE IF RETREATING WOULD BE NEEDED WITH A
C     LINEAR MODEL.  ONLY THEN RETREAT.
      IF (RC.LE.ONE .OR. .NOT. USEQ) GO TO 750
C     DO(EVALUATE LINEAR MODEL)
      NV = NVARS
      ASSIGN 740 TO IGOELM
      GO TO 940

  740 CONTINUE
      PVL = DNRM2(MIN(MEQUA,NVARS+1),WJ(MCONST+1,NV+1),1)
      IF (PVL.GT.PB) RC = 4.D0* (PVL-PB)/ (PVL+PD)
  750 CONTINUE
      RG = MAX(RG,RC)
      IF (IPRINT.GT.0) THEN
          WRITE (NOUT,9011) ITERS,FC,PV,RC,AJN,K,KL,FB,ALPHA,BBOOST,
     .      NIT,USEQ,NTTERM
          WRITE (NOUT,9001) '  X=', (X(J),J=1,NVARS)
          WRITE (NOUT,9001) ' DX=', (DX(J),J=1,NVARS)
          WRITE (NOUT,9001) '  B=', (B(J),J=1,NVARS)
          WRITE (NOUT,9001) ' LB=', (BLB(J),J=1,NALL)
          WRITE (NOUT,9001) ' UB=', (BUB(J),J=1,NALL)
          WRITE (NOUT,'('' ///END OF ITERATION.///'')')
      END IF

      RETREA = RC .GT. 1
      IF ( .NOT. RETREA) THEN
          CHG = ONE
          T2 = ZERO
          DO 810 J = 1,NVARS
             BOLD = B(J)
             T = DX(J)/BOLD
             ALB = TWO
             AUB = TWO
C    IF USER GIVES BOUNDS, AND THESE BOUNDS ARE HIT,
C    DO NOT DETRACT FROM DECLARING A FULL NEWTON STEP.
             ICASE = IND(J)
             GO TO (760,770,780,790),ICASE
C     CASE 1
  760        ALB = (X(J)-BL(J))/BOLD
             AUB = -SEMIBG
             GO TO 800
C     CASE 2
  770        AUB = (X(J)-BU(J))/BOLD
             ALB = -SEMIBG
             GO TO 800
C     CASE 3
  780        ALB = (X(J)-BL(J))/BOLD
             AUB = (X(J)-BU(J))/BOLD
             GO TO 800
C     CASE 4
  790        CONTINUE
             ALB = -SEMIBG
             AUB = -SEMIBG
C     END CASE
  800        CONTINUE
             IF (T.EQ.ONE) THEN
                 T2 = ONE
                 B(J) = BOLD + BOLD
                 CHG = CHG*CHGFAC

             ELSE
                 IF (ABS(T).LT..25D0 .AND. DX(J).NE.ZERO) THEN
                     B(J) = SIGN(.25D0*BOLD,DX(J)) + 3.D0*DX(J)

                 ELSE
                     B(J) = SIGN(BOLD,DX(J))
                 END IF

             END IF
C     THIS TEST AVOIDS THE USER BOUNDS IN DECLARING A NEWTON STEP.
             IF (ABS(ALB-T).GE..01D0*ABS(T) .AND. ABS(AUB-T).GE..01D0*
     .           ABS(T)) THEN
                 IF (T.GT.ZERO) THEN
                     T2 = MAX(T2,T)

                 ELSE
                     T2 = MAX(T2,-T/C1516)
                 END IF

             END IF

  810     CONTINUE
          FULNWT = T2 .LT. .99D0
          DXNRM = ABS(DX(IDAMAX(NVARS,DX,1)))
C     DO BLOCK
C     TEST FOR SMALL ABSOLUTE CHANGE IN X VALUES.
          TERM = DXNRM .LT. TOLD .AND. FULNWT
          IF (TERM) THEN
              IGO = 6
C     VALUE MEANS CHANGE (IN PARAMETERS) WAS SMALL AND A
C     FULL STEP (NOT HITTING TRUST CONSTRAINTS) TAKEN.
C     EXIT BLOCK
              GO TO 820

          END IF

          TERM = DXNRM .LT. DNRM2(NVARS,X,1)*TOLX .AND. FULNWT
          TERM = TERM .AND. (ITERS.GT.1)
          IF (TERM) THEN
              IGO = 7
C     VALUE MEANS RELATIVE CHANGE IN PARAMETERS WAS SMALL AND A
C     FULL STEP (NOT HITTING CONSTRAINTS WITH AT LEAST 2 ITERATIONS)
C     WAS TAKEN.
C     EXIT BLOCK
              GO TO 820

          END IF
C     EXIT IF
          GO TO 830
C     END BLOCK
  820     CONTINUE
C     CYCLE FOREVER
          GO TO 30

      END IF

  830 CONTINUE
      FL = FC
C     END FOREVER
      GO TO 30

  840 CONTINUE
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 1020
C     PROCEDURE(EVALUATE QUADRATIC MODEL)
  850 CONTINUE
C     IF THE MODEL IS GENUINELY QUADRATIC, ADD IN THE EXTRA
C     TERMS AND COMPUTE THE SECOND DERIVATIVE INFORMATION.
      IF (USEQ .AND. .NOT. LINMOD) THEN
C     COMPUTE THE DOT PRODUCT OF CURRENT PROPOSED STEP AND
C     PAST DIRECTIONS REPRESENTED IN THE MODEL.
          DO 870 L = 1,NP - 1
             T = ZERO
             DO 860 J = 1,NVARS
                T = T + DX(J)* (XP(J,L+1)-XP(J,1))
  860        CONTINUE
             PJ(L) = T
  870     CONTINUE
C     STORAGE LAYOUT, WITH K = J**T, OF WJ(*,*).
C     [J    :  I    : F+R ]
C     [H    :  K    : K*R ]
C     ADD IN THE QUADRATIC TERMS FOR THE FUNCTION.
          DO 880 L = 1,NP - 1
             JK = MIN(NVARS+L+1,MEQUA)
             CALL DAXPY(JK,.5D0*PJ(L)**2,QC(1,L+1),1,WJ(MCONST+1,NV+1),
     .                  1)
  880     CONTINUE
C     ADD THE LINEAR TERMS TO THE INNER LOOP JACOBIAN.
          DO 900 L = 1,NP - 1
             JK = MIN(NVARS+L+1,MEQUA)
             DO 890 J = 1,NVARS
                CALL DAXPY(JK,PJ(L)* (XP(J,L+1)-XP(J,1)),QC(1,L+1),1,
     .                     WJ(MCONST+1,J),1)
  890        CONTINUE
  900     CONTINUE
C     COMPUTE THE UPPER TRIANGULAR PART OF THE SECOND
C     DERIVATIVE TERMS.
          DO 930 I = 1,NVARS
             DO 920 J = I,NVARS
                DO 910 L = 1,NP - 1
                   JK = MIN(NVARS+L+1,MEQUA)
                   WJ(MCONST+MK+I,J) = WJ(MCONST+MK+I,J) +
     .                                 (XP(J,L+1)-XP(J,1))*
     .                                 (XP(I,L+1)-XP(I,1))*
     .                                 DDOT(JK,DX(NVARS+1),1,QC(1,L+1),
     .                                 1)
  910           CONTINUE
  920        CONTINUE
  930     CONTINUE
      END IF
C     END PROCEDURE
      GO TO IGOEQM
C     PROCEDURE(EVALUATE LINEAR MODEL)
  940 CONTINUE
C     TRANSFER THE JACOBIAN THAT WOULD RESULT FROM
C     USING JUST A LINEAR MODEL.
      DO 950 J = 1,NVARS
         IF (JACTRI) THEN
             JK = J

         ELSE
             JK = MEQUA
         END IF

         CALL DCOPY(JK,FJAC(MCON+1,J),1,WJ(MCONST+1,J),1)
  950 CONTINUE
C     TRANSFER THE PRESENT VALUES OF THE FUNCTION.
      CALL DCOPY(MIN(MEQUA,NVARS+1),FJAC(MCON+1,NVARS+1),1,
     .           WJ(MCONST+1,NV+1),1)
C     CHANGE SIGN FOR THE MODEL PROBLEM.
      DO 960 I = 1,MIN(MEQUA,NVARS+1)
         WJ(MCONST+I,NV+1) = -WJ(MCONST+I,NV+1)
  960 CONTINUE
C     COMPUTE THE LINEAR TERM OF THE MODEL.
      DO 970 J = 1,NVARS
         IF (JACTRI) THEN
             JK = J

         ELSE
             JK = MEQUA
         END IF

         CALL DAXPY(JK,DX(J),WJ(MCONST+1,J),1,WJ(MCONST+1,NV+1),1)
  970 CONTINUE
C     END PROCEDURE
      GO TO IGOELM
C     PROCEDURE(TEST FOR CONVERGENCE)
  980 CONTINUE
      TERM = ITERS .GE. ITMAX
      IF (TERM) THEN
          IGO = 8
C     VALUE MEANS THAT MAX. NUMBER OF ALLOWED ITERATIONS TAKEN.
C     EXIT PROCEDURE
          GO TO 1000

      END IF
C     TEST FOR SMALL FUNCTION NORM.
      TERM = FC .LE. TOLF
C     IF HAVE CONSTRAINTS MUST ALLOW AT LEAST ONE MOVE.
      TERM = TERM .AND. (MCON.EQ.0 .OR. ITERS.GT.1)
      IF (TERM) THEN
          IGO = 2
C     VALUE MEANS FUNCTION NORM WAS SMALL.
C     EXIT PROCEDURE
          GO TO 1000

      END IF
C     DO(TEST FOR NO CHANGE)
      GO TO 1010

  990 CONTINUE
      TERM = TERM .AND. .NOT. RETREA
      IF (TERM) THEN
          IGO = 3
C     VALUE MEANS THE FUNCTION IS PROBABLY REACHING A LOCAL MINIMUM
C     BUT MOVES ARE STILL HITTING TRUST REGION CONSTRAINTS.
          IF (FULNWT) IGO = 4
C     VALUE MEANS THAT FUNCTION IS REACHING A LOCAL MINIMUM
C     AND MOVES ARE NOT HITTING THE TRUST REGION CONSTRAINTS.
          IF (IGO.EQ.3) TERM = TERM .AND. .NOT. MUSTCN
C     EXIT PROCEDURE
          GO TO 1000

      END IF
C     END PROCEDURE
 1000 CONTINUE
      GO TO IGOTFC
C     PROCEDURE(TEST FOR NO CHANGE)
 1010 CONTINUE
      TERM = (ABS(FB-PV).LE.TOLSNR*FB) .AND. (ABS(FC-PV).LE.FB*TOLP)
      TERM = TERM .AND. (ABS(FC-FL).LE.FB*TOLSNR)
      TERM = TERM .AND. (ABS(PVL-PV).LE.FB*TOLSNR)
C     END PROCEDURE
      GO TO 990
C     PROCEDURE(RETURN TO USER PROGRAM UNIT)
 1020 CONTINUE
      RETURN
C     PROCEDURE(INITIALIZE OTHER VALUES)
 1030 CONTINUE
C     THE NUMBER OF PAST DIFFERENCES USED IN THE QUADRATIC MODEL.
      NP = 0
C     IF NO MORE EQUATIONS THAN VARIABLES, NO NEED TO
C     PRETRIANGULARIZE THE JACOBIAN MATRIX.
      JACTRI = MEQUA .GT. NVARS
C     MAKE SURE THAT VARIABLES SATISFY CONSTRAINTS.
C     GENERALLY THIS MAY TAKE A CALL TO DBOCLS().
C     AS LONG AS THE FUNCTIONS ARE DEFINED AT POINTS
C     THAT DO NOT SATISFY THE CONSTRAINTS, THE FIRST
C     ALGORITHM STEP WILL BRING IT ONTO THE CONSTRAINTS.
C
      DO 1090 J = 1,NVARS
         GO TO (1040,1050,1060,1070),IND(J)
C     CASE 1
 1040    X(J) = MAX(X(J),BL(J))
         GO TO 1080
C     CASE 2
 1050    X(J) = MIN(X(J),BU(J))
         GO TO 1080
C     CASE 3
 1060    X(J) = MAX(X(J),BL(J))
         X(J) = MIN(X(J),BU(J))
         GO TO 1080
C     CASE 4
 1070    CONTINUE
C     END CASE
 1080    CONTINUE
 1090 CONTINUE
      ITERS = 0
      NALL = MCON + NVARS
      CHGFAC = TWO** (-ONE/REAL(NVARS))
      C1516 = 15.D0/16.D0
      SEMIBG = 1.D10
C     END PROCEDURE
      GO TO 20
C     PROCEDURE(PROCESS OPTION ARRAY)
 1100 CONTINUE
      IPRINT = 0
C     I1MACH(2)=STANDARD OUTPUT UNIT.
      NOUT = I1MACH(2)
C     D1MACH(4)=RELPR=MACHINE REL. PREC.
      T = D1MACH(4)
      TOLF = SQRT(T)
      TOLUSE = TOLF
      TOLX = TOLF
      TOLD = TOLF
      TOLSNR = 1.D-5
      TOLP = 1.D-5
      COND = 30.
      ITMAX = 75
      LEVEL = 1
      IPLS = 0
      PASSB = .FALSE.
      NOQUAD = .FALSE.
      REVERS = .FALSE.
      MUSTCN = .FALSE.
      LP = 1
      LPDIFF = 0
C     DO FOREVER
 1110 CONTINUE
      LP = LP + LPDIFF
      LPDIFF = 2
      KP = IOPT(LP)
      NEWOPT = KP .GT. 0
      JP = ABS(KP)
C     SEE IF THIS IS THE LAST OPTION..
      IF (JP.EQ.99) THEN
          IF (NEWOPT) THEN
C     THE POINTER TO THE START OF OPTIONS FOR THE INNER LOOP
C     SOLVER MUST SATISFY THE REQUIREMENTS FOR THAT OPTION ARRAY.
              IF (IPLS.EQ.0) IPLS = LP
C     EXIT FOREVER
              GO TO 1120
*
          ELSE
              LPDIFF = 1
C     CYCLE FOREVER
              GO TO 1110

          END IF

      END IF
C     CHANGE PRINT OPTION.
      IF (JP.EQ.1) THEN
          IF (NEWOPT) IPRINT = IOPT(LP+1)
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     SEE IF MAX. NUMBER OF ITERATIONS CHANGING.
      IF (JP.EQ.2) THEN
          IF (NEWOPT) ITMAX = IOPT(LP+1)
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     SEE IF BOUNDS FOR THE TRUST REGION ARE BEING PASSED.
      IF (JP.EQ.3) THEN
          IF (NEWOPT) THEN
              CALL DCOPY(NVARS,ROPT(IOPT(LP+1)),1,BB,1)
              PASSB = .TRUE.
          END IF
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     CHANGE TOLERANCE ON THE LENGTH OF THE RESIDUALS.
      IF (JP.EQ.4) THEN
          IF (NEWOPT) TOLF = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     CHANGE TOLERANCE ON THE NORM OF THE RELATIVE
C     CHANGE TO THE PARAMETERS.
      IF (JP.EQ.5) THEN
          IF (NEWOPT) TOLX = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     CHANGE TOLERANCE ON ABSOLUTE CHANGE TO THE PARAMETERS.
      IF (JP.EQ.6) THEN
          IF (NEWOPT) TOLD = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF

      IF (JP.EQ.7) THEN
C     CHANGE TOLERANCE FOR RELATIVE AGREEMENT BETWEEN
C     BEST FUNCTION NORM, LAST FUNCTION NORM AND THE
C     CURRENT FUNCTION NORM.
          IF (NEWOPT) TOLSNR = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF

      IF (JP.EQ.8) THEN
C     CHANGE TOLERANCE FOR AGREEMENT BETWEEN PREDICTED
C     VALUE OF RESIDUAL NORM AND THE PREVIOUS VALUE OF
C     THE RESIDUAL NORM.
          IF (NEWOPT) TOLP = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF

      IF (JP.EQ.9) THEN
C     CHANGE TOLERANCE SUCH THAT RELATIVE CHANGES IN THE
C     VALUES OF THE PARAMETERS IMPLY THAT THE PREVIOUS
C     VALUE OF THE FUNCTION WILL NOT BE USED IN THE
C     QUADRATIC MODEL.
          IF (NEWOPT) TOLUSE = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF

      IF (JP.EQ.10) THEN
C     CHANGE THE LARGEST CONDITION NUMBER TO ALLOW WHEN
C     SOLVING FOR THE QUADRATIC COEFFICIENTS OF THE MODEL.
          IF (NEWOPT) COND = ROPT(IOPT(LP+1))
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     CHANGE THE PRINT LEVEL IN THE ERROR PROCESSOR.
      IF (JP.EQ.11) THEN
          IF (NEWOPT) LEVEL = IOPT(LP+1)
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     PASS AN OPTION ARRAY TO THE CONSTRAINED LINEAR SOLVER.
C     THIS OPTION IS A POINTER TO THE START OF THE OPTION
C     ARRAY FOR THE SUBPROGRAM.
      IF (JP.EQ.12) THEN
          IF (NEWOPT) IPLS = IOPT(LP+1)
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     MOVE THE PROCESSING POINTER BY THE VALUE IN THE
C     NEXT ENTRY OF THE OPTION ARRAY.  THIS DEVICE IS
C     INCLUDED SO THAT PASSING OPTIONS TO LOWER LEVEL
C     SUBROUTINES IS EASY TO DO.
      IF (JP.EQ.13) THEN
          IF (NEWOPT) LPDIFF = IOPT(LP+1)
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     OPTION TO SUPPRESS USING THE QUADRATIC MODEL, EVER.
      IF (JP.EQ.14) THEN
          IF (NEWOPT) NOQUAD = IOPT(LP+1) .EQ. 1
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     MORE STORAGE WAS GIVEN FOR THE QUADRATIC MODEL ARRAYS.
C     THIS OPTION WAS PROCESSED BY THE INTERFACE UNIT.
C     IF(JP.EQ.15)CYCLE FOREVER
      IF (JP.EQ.15) GO TO 1110
C     USE FORWARD COMMUNICATION TO GET THE DERIVATIVES
C     AND FUNCTION VALUES.
      IF (JP.EQ.16) THEN
          IF (NEWOPT) REVERS = IOPT(LP+1) .EQ. 1
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     FORCE A FULL NEWTON STEP WHEN NEAR THE MINIMUM.
C     DO NOT ALLOW CONVERGENCE CLAIMS WHEN HITTING BOUNDS.
      IF (JP.EQ.17) THEN
          IF (NEWOPT) MUSTCN = IOPT(LP+1) .EQ. 1
C     CYCLE FOREVER
          GO TO 1110

      END IF
C     SAW AN OPTION (OR GARBAGE) THAT IS NOT ON THE LIST.
      XMESS =
     .'DQEDMN. INVALID OPTION PROCESSED. I1=IOPT(*) ENTRY. I2=IOPT(I1).'
      NERR = 07
      IGO = 15
      CALL CHRCNT(XMESS,NMESS)
      CALL XERRWV(XMESS,NMESS,NERR,LEVEL,2,LP,IOPT(LP),0,RDUM,RDUM)
      IFLAG = 0
C     DO(RETURN TO USER PROGRAM UNIT)
      GO TO 1020
C     END FOREVER
C     END PROCEDURE
 1120 CONTINUE
      GO TO 10

 9001 FORMAT (A4,1P10D12.4/ (4X,10D12.4))
 9011 FORMAT ('0ITER.=',I3,' FC=',1PD10.4,' PV=',1PD10.4,03X,' RC=',
     .  1PD10.4,' J**T*F=',1PD10.4,/,' K=',I4,' KL=',I4,/10X,' FB=',
     .  1PD10.4,' AL=',1PD10.4,' BB=',1PD12.4/' INNER ITERATIONS =',I5,
     .  ' USE QUAD. MODEL?=',L5,' NUM. OF TERMS =',I5)
 9021 FORMAT (' MODEL RESIDUAL.GE.CURRENT F. QUITTING.',1P2D12.5)
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
      do i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
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
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision SA, the multiplier.
c
c    Input/output, double precision X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      double precision da,dx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      end do
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      end do
      if( n .lt. 5 ) return
   40 continue
      do i = m+1, n, 5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      end do

      return
      end
      subroutine dswap ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DSWAP interchanges two vectors.
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
c    Input/output, double precision X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, double precision Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n
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
      do i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
      end do

      return
      end
      SUBROUTINE DVOUT(N,DX,IFMT,IDIGIT)
C     REVISED FEB. 27, 1981.
      DOUBLE PRECISION DX(*)
      CHARACTER IFMT*(*)
C
C     DOUBLE PRECISION VECTOR OUTPUT ROUTINE.
C
C  INPUT..
C
C  N,DX(*) PRINT THE DOUBLE PRECISION ARRAY DX(I),I=1,...,N, ON
C          OUTPUT UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT
C          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
C          STEP. THE COMPONENTS DX(I) ARE INDEXED, ON OUTPUT,
C          IN A PLEASANT FORMAT.
C  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
C          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT
C                WRITE(LOUT,IFMT)
C  IDIGIT  PRINT AT LEAST IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
C          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 6,14,20 OR 28
C          WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
C          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
C          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY DX(*). (THIS
C          CAN BE USED ON MOST TIME-SHARING TERMINALS). IF
C          IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN
C          BE USED ON MOST LINE PRINTERS).
C
C  EXAMPLE..
C
C  PRINT AN ARRAY CALLED (COSTS OF PURCHASES) OF LENGTH 100 SHOWING
C  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
C  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
C
C     DOUBLE PRECISION COSTS(100)
C     N = 100
C     IDIGIT = -6
C     CALL DVOUT(N,COSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
C
C
C
C     AUTHORS    JOHN A. WISNIEWSKI   SANDIA LABS ALBUQUERQUE.
C                RICHARD J. HANSON    SANDIA LABS ALBUQUERQUE.
C     DATE       JULY 27,1978.
C
C
C     GET THE UNIT NUMBER WHERE OUTPUTWILL BE WRITTEN.
      J=2
      LOUT=I1MACH(J)
      WRITE(LOUT,IFMT)
      IF(N.LE.0) RETURN
      NDIGIT = IDIGIT
      IF(IDIGIT.EQ.0) NDIGIT = 6
      IF(IDIGIT.GE.0) GO TO 80
C
      NDIGIT = -IDIGIT
      IF(NDIGIT.GT.6) GO TO 20
C
      DO 10 K1=1,N,4
      K2 = MIN0(N,K1+3)
      WRITE(LOUT,1000) K1,K2,(DX(I),I=K1,K2)
   10 CONTINUE
      RETURN
C
   20 CONTINUE
      IF(NDIGIT.GT.14) GO TO 40
C
      DO 30 K1=1,N,2
      K2 = MIN0(N,K1+1)
      WRITE(LOUT,1001) K1,K2,(DX(I),I=K1,K2)
   30 CONTINUE
      RETURN
C
   40 CONTINUE
      IF(NDIGIT.GT.20) GO TO 60
C
      DO 50 K1=1,N,2
      K2=MIN0(N,K1+1)
      WRITE(LOUT,1002) K1,K2,(DX(I),I=K1,K2)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      DO 70 K1=1,N
      K2 = K1
      WRITE(LOUT,1003) K1,K2,(DX(I),I=K1,K2)
   70 CONTINUE
      RETURN
C
   80 CONTINUE
      IF(NDIGIT.GT.6) GO TO 100
C
      DO 90 K1=1,N,8
      K2 = MIN0(N,K1+7)
      WRITE(LOUT,1000) K1,K2,(DX(I),I=K1,K2)
   90 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(NDIGIT.GT.14) GO TO 120
C
      DO 110 K1=1,N,5
      K2 = MIN0(N,K1+4)
      WRITE(LOUT,1001) K1,K2,(DX(I),I=K1,K2)
  110 CONTINUE
      RETURN
C
  120 CONTINUE
      IF(NDIGIT.GT.20) GO TO 140
C
      DO 130 K1=1,N,4
      K2 = MIN0(N,K1+3)
      WRITE(LOUT,1002) K1,K2,(DX(I),I=K1,K2)
  130 CONTINUE
      RETURN
C
  140 CONTINUE
      DO 150 K1=1,N,3
      K2 = MIN0(N,K1+2)
      WRITE(LOUT,1003) K1,K2,(DX(I),I=K1,K2)
  150 CONTINUE
      RETURN
 1000 FORMAT(1X,I4,' - ',I4,1X,1P8D14.5)
 1001 FORMAT(1X,I4,' - ',I4,1X,1P5D22.13)
 1002 FORMAT(1X,I4,' - ',I4,1X,1P4D28.19)
 1003 FORMAT(1X,I4,' - ',I4,1X,1P3D36.27)
      END
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
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer
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
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the index of element having maximum absolute value.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer IDAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      double precision dx(*),dmax
      integer idamax
      integer i,incx,ix,n

      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do  i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
      end do
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do i = 2,n
        if( dmax .lt. dabs(dx(i)) ) then
          idamax = i
          dmax = dabs(dx(i))
        end if
      end do

      return
      end
      SUBROUTINE IVOUT(N,IX,IFMT,IDIGIT)
C     REVISED FEB. 27, 1981.
      DIMENSION IX(N)
      CHARACTER IFMT*(*)
C
C     INTEGER VECTOR OUTPUT ROUTINE.
C
C  INPUT..
C
C  N,IX(*) PRINT THE INTEGER ARRAY IX(I),I=1,...,N, ON OUTPUT
C          UNIT LOUT. THE HEADING IN THE FORTRAN FORMAT
C          STATEMENT IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST
C          STEP. THE COMPONENTS IX(I) ARE INDEXED, ON OUTPUT,
C          IN A PLEASANT FORMAT.
C  IFMT(*) A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON OUTPUT
C          UNIT LOUT WITH THE VARIABLE FORMAT FORTRAN STATEMENT
C                WRITE(LOUT,IFMT)
C  IDIGIT  PRINT UP TO IABS(IDIGIT) DECIMAL DIGITS PER NUMBER.
C          THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,10 OR 14
C          WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF
C          PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE UTILIZED
C          TO WRITE EACH LINE OF OUTPUT OF THE ARRAY IX(*). (THIS
C          CAN BE USED ON MOST TIME-SHARING TERMINALS). IF
C          IDIGIT.GE.0, 133 PRINTING COLUMNS ARE UTILIZED. (THIS CAN
C          BE USED ON MOST LINE PRINTERS).
C
C  EXAMPLE..
C
C  PRINT AN ARRAY CALLED (COSTS OF PURCHASES) OF LENGTH 100 SHOWING
C  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING
C  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE.
C
C     DIMENSION ICOSTS(100)
C     N = 100
C     IDIGIT = -6
C     CALL IVOUT(N,ICOSTS,'(''1COSTS OF PURCHASES'')',IDIGIT)
C
C
C
C     AUTHORS    JOHN A. WISNIEWSKI   SANDIA LABS ALBUQUERQUE.
C                RICHARD J. HANSON    SANDIA LABS ALBUQUERQUE.
C     DATE       JULY 27,1978.
C
C
C     GET THE UNIT NUMBER WHERE OUTPUTWILL BE WRITTEN.
      J=2
      LOUT=I1MACH(J)
      WRITE(LOUT,IFMT)
      IF(N.LE.0) RETURN
      NDIGIT = IDIGIT
      IF(IDIGIT.EQ.0) NDIGIT = 4
      IF(IDIGIT.GE.0) GO TO 80
C
      NDIGIT = -IDIGIT
      IF(NDIGIT.GT.4) GO TO 20
C
      DO 10 K1=1,N,10
      K2 = MIN(N,K1+9)
      WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   10 CONTINUE
      RETURN
C
   20 CONTINUE
      IF(NDIGIT.GT.6) GO TO 40
C
      DO 30 K1=1,N,7
      K2 = MIN(N,K1+6)
      WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
   30 CONTINUE
      RETURN
C
   40 CONTINUE
      IF(NDIGIT.GT.10) GO TO 60
C
      DO 50 K1=1,N,5
      K2=MIN(N,K1+4)
      WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
   50 CONTINUE
      RETURN
C
   60 CONTINUE
      DO 70 K1=1,N,3
      K2 = MIN(N,K1+2)
      WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
   70 CONTINUE
      RETURN
C
   80 CONTINUE
      IF(NDIGIT.GT.4) GO TO 100
C
      DO 90 K1=1,N,20
      K2 = MIN(N,K1+19)
      WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   90 CONTINUE
      RETURN
C
  100 CONTINUE
      IF(NDIGIT.GT.6) GO TO 120
C
      DO 110 K1=1,N,15
      K2 = MIN(N,K1+14)
      WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
  110 CONTINUE
      RETURN
C
  120 CONTINUE
      IF(NDIGIT.GT.10) GO TO 140
C
      DO 130 K1=1,N,10
      K2 = MIN(N,K1+9)
      WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
  130 CONTINUE
      RETURN
C
  140 CONTINUE
      DO 150 K1=1,N,7
      K2 = MIN(N,K1+6)
      WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
  150 CONTINUE
      RETURN
 1000 FORMAT(1X,I4,' - ',I4,20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,7(1X,I15))
      END
      SUBROUTINE XERROR (XMESS, NMESS, NERR, LEVEL)
C
C     REPLACEMENT SUBROUTINE FOR THE SLATEC XERROR PACKAGE SUBROUTINE.
      CHARACTER*120 XMESS
      IF (LEVEL.GE.1) THEN
      IERR=I1MACH(3)
      WRITE(IERR,'(1X,A)') XMESS(1:NMESS)
      WRITE(IERR,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)')
     .      NERR,LEVEL
      END IF
      RETURN
      END
      SUBROUTINE XERRWV(XMESS,NMESS,NERR,LEVEL,NI,I1,I2,
     .           NR,R1,R2)
C
C     REPLACEMENT SUBROUTINE FOR THE SLATEC XERRWV PACKAGE SUBROUTINE.
      CHARACTER*120 XMESS
      IF(LEVEL.LT.1) RETURN
      IERR=I1MACH(3)
      WRITE(IERR,'(1X,A)') XMESS(1:NMESS)
      WRITE(IERR,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)')
     .      NERR,LEVEL
      ASSIGN 80 TO IGOIPR
      GO TO (10,20),NI
      GO TO 80
   10 CONTINUE
C
C     WRITE ONE INTEGER VALUE.
      WRITE(IERR,'('' I1 = '',I8)') I1
      GO TO IGOIPR
   20 ASSIGN 30 TO IGOIPR
      GO TO 10
   30 CONTINUE
C
C     WRITE ONE INTEGER VALUE.
      WRITE(IERR,'('' I2 = '',I8)') I2
   80 CONTINUE
      ASSIGN 40 TO IGOIPR
      GO TO (50,60),NR
   40 RETURN
   50 CONTINUE
C
C     WRITE ONE REAL VALUE.
      WRITE(IERR,'('' R1 = '',1PE15.7)') R1
      GO TO IGOIPR
   60 ASSIGN 70 TO IGOIPR
      GO TO 50
   70 CONTINUE
C
C     WRITE REAL VALUE.
      WRITE(IERR,'('' R2 = '',1PE15.7)') R2
      GO TO 40
      END
