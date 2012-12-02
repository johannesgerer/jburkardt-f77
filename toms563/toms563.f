      SUBROUTINE CL1(NEQNS,NEQC,NIQC,NVARS,NACT,IFL,MXS,PSW,
     *               E,NER,X,F,EL1N,RES,INDX,W)
C
      INTEGER IFL,INDX(1),MXS,NACT,NEQC,NEQNS,NER,NIQC,NVARS
      LOGICAL PSW
      REAL E(NER,1),EL1N,F(1),RES(1),W(1),X(1)
C
C
C     **********************************************************
C     A PROGRAM FOR THE SOLUTION IN THE   L1  SENSE
C     OF A LINEAR-EQUATION SYSTEM
C     (WITH OR WITHOUT LINEAR CONSTRAINTS).
C     RICHARD H. BARTELS AND ANDREW R. CONN.
C     LATEST UPDATE .... 13 APRIL, 1980.
C
C     DEVELOPMENT OF THIS PROGRAM WAS SUPPORTED
C     IN PART BY U. S. NSF GRANT DCR75-07817
C     AND BY FUNDS FROM THE NATIONAL BUREAU OF STANDARDS,
C     AND IN PART BY CANADIAN NRC GRANT A8639.
C
C     +++++ PARAMETERS +++++
C     ----------------------------------------------------------
C                           INPUT
C     NAME  TYPE  SUBSCRPT  OUTPUT        DESCRIPTION
C                           SCRATCH
C     ..........................................................
C     NEQNS INT.    NONE      IN      NUMBER OF EQUATIONS
C                                     (MAY BE ZERO)
C
C     NEQC  INT.    NONE      IN      NUMBER OF EQUALITY
C                                     CONSTRAINTS
C                                     (MAY BE ZERO)
C
C     NIQC  INT.    NONE      IN      NUMBER OF INEQUALITY
C                                     CONSTRAINTS
C                                     (MAY BE ZERO)
C
C     NVARS INT.    NONE      IN      NUMBER OF VARIABLES
C
C     NACT  INT.    NONE      OUT     NUMBER OF ACTIVE
C                                     EQUATIONS/CONSTRAINTS
C                                     AT TERMINATION
C                                     (IF ANY, THEIR ASSOCIATED
C                                      COLUMN POSITIONS IN  E  WILL
C                                      BE LISTED IN  INDX(1)
C                                      THROUGH  INDX(NACT) )
C
C     IFL   INT.    NONE      OUT     TERMINATION CODE
C                                     (SEE BELOW)
C
C     MXS   INT.    NONE      IN      MAXIMUM NUMBER OF STEPS
C                                     ALLOWED
C
C     PSW   LOGIC.  NONE      IN      PRINT SWITCH
C                                     (SEE BELOW)
C
C     E     REAL     2        IN      EQUATION/CONSTRAINT MATRIX
C                                     THE FIRST  NEQNS  COLUMNS
C                                     (SEE NOTE BELOW) SPECIFY
C                                     EQUATIONS, THE REMAINING
C                                     COLUMNS (IF ANY) SPECIFY
C                                     CONSTRAINTS.
C
C     NER   INT.    NONE      IN      ROW DIMENSION OF E
C
C     X     REAL     1        IN      STARTING VALUES FOR THE
C                                     UNKNOWNS (USE ZEROS IF NO
C                                     GUESS IS AVAILABLE)
C                             OUT     TERMINATION VALUES FOR
C                                     THE UNKNOWNS
C
C     F     REAL     1        IN      EQUATION/CONSTRAINT
C                                     RIGHT-HAND SIDES
C
C     EL1N  REAL    NONE      OUT     L1 NORM OF EQUATION
C                                     RESIDUALS AT TERMINATION
C
C     RES   REAL     1        OUT     EQUATION/CONSTRAINT
C                                     RESIDUALS AT TERMINATION
C
C     INDX  INT.     1        OUT     INDEX VECTOR USED TO RECORD
C                                     THE ORDER IN WHICH THE COLUMNS
C                                     OF  E  ARE BEING PROCESSED
C
C     W     REAL     1        SCR.    WORKING STORAGE
C     ----------------------------------------------------------
C
C     +++++ PURPOSE +++++
C     ----------------------------------------------------------
C     THIS SUBROUTINE SOLVES THE   NEQNS BY NVARS
C     SYSTEM OF EQUATIONS
C
C                       (A-TRANSPOSE) * X   ==   B
C
C     SUBJECT TO THE  NEQC   CONSTRAINTS
C
C                       (G-TRANSPOSE) * X  .EQ.  H
C
C     AND THE  NIQC  INEQUALITY CONSTRAINTS
C
C                       (C-TRANSPOSE) * X  .GE.  D
C
C     FOR THE UNKNOWNS  X(1),...,X(NVARS).
C
C     THE PROBLEM MUST BE WELL-POSED, NONTRIVIAL
C     AND OVERDETERMINED IN THE SENSE THAT
C
C                          NVARS .GE. 1
C                          NEQNS .GE. 0
C                          NEQC  .GE. 0
C                          NIQC  .GE. 0
C               NEQNS+NEQC+NIQC  .GE. NVARS.
C
C     FURTHER, NO COLUMN OF  A, G  OR  C  SHOULD BE ZERO.
C     IF THESE CONDITIONS ARE NOT MET, THE PROGRAM
C     WILL TERMINATE WITHOUT PERFORMING ANY SUBSTANTIVE
C     COMPUTATIONS.
C
C     A POINT  X  IS A SOLUTION IF IT MINIMIZES THE EQUATION
C     RESIDUALS FROM AMONG ALL POINTS WHICH SATISFY THE
C     CONSTRAINTS.  AT ANY (NONDEGENERATE) SOLUTION
C     THERE WILL BE  NACT  EQUATIONS AND CONSTRAINTS
C     WHOSE RESIDUALS
C
C          (A(I)-TRANSPOSE) * X - B(I)
C
C          (G(I)-TRANSPOSE) * X - H(I)
C
C     AND
C
C          (C(I)-TRANSPOSE) * X - D(I)
C
C     ARE ZERO.
C
C     THE COLUMNS OF  (A,G,C)  CORRESPONDING TO THE ZERO RESIDUALS
C     ARE REFERRED TO AS  ACTIVE COLUMNS  THROUGHOUT THIS LISTING.
C     THE NUMBERS OF THE ACTIVE COLUMNS ARE MAINTAINED AS THE
C     ENTRIES  1,...,NACT  OF THE ARRAY  INDX.
C
C     A SOLUTION  X  IS FOUND BY MINIMIZING A PIECEWISE
C     LINEAR PENALTY FUNCTION FORMED FROM THE  L1
C     NORM OF THE EQUATION RESIDUALS AND THE SUM OF THE
C     INFEASIBILITIES IN THE CONSTRAINTS.
C     THE MINIMIZATION PROCEEDS IN A STEP-BY-STEP
C     FASHION, TERMINATING AFTER A FINITE NUMBER OF STEPS.
C
C     NOTE THAT  A, G  AND  C  APPEAR TRANSPOSED IN THE
C     PROBLEM FORMULATION.  HENCE IT IS THE COLUMNS OF  (A,G,C)
C     WHICH DEFINE THE EQUATIONS AND CONSTRAINTS RESPECTIVELY.
C
C     THE ARRAY  E  IS A COMPOSITE OF   A, G AND C
C     AND  F  IS A COMPOSITE OF  B, H  AND  D.
C     E  SHOULD CONTAIN  A  AS ITS FIRST  NEQNS  COLUMNS.
C     IT SHOULD CONTAIN  G  AS ITS NEXT  NEQC  COLUMNS AND
C     CONTAIN  C  AS ITS REMAINING  NIQC  COLUMNS.
C     SIMILARLY  F  SHOULD CONTAIN  B  AS ITS FIRST
C     NEQNS  COMPONENTS,  H  AS ITS NEXT  NEQC  COMPONENTS
C     AND  D  AS ITS LAST  NIQC  COMPONENTS.
C     ----------------------------------------------------------
C
C     +++++ ARRAYS +++++
C     ----------------------------------------------------------
C     E  IS TO BE DIMENSIONED AT LEAST    N  BY  M,
C     X                       AT LEAST    N,
C     F                       AT LEAST    M,
C     RES                     AT LEAST    M,
C     INDX                    AT LEAST    M,
C     W                       AT LEAST    ((3*N*N+11*N+2)/2) + (2*M).
C
C                                         WHERE  N = NVARS  AND
C                                         M = NEQNS+NEQC+NIQC
C     ----------------------------------------------------------
C
C     +++++ INITIALIZATION +++++
C     ----------------------------------------------------------
C     THE USER MUST INITIALIZE
C
C          NEQNS,NEQC,NIQC,NVARS,MXS,PSW,E,NER,X,F .
C
C     THE FOLLOWING ARE SET BY  CL1
C     AND DO NOT REQUIRE INITIALIZATION
C
C          NACT,INDX,RES .
C
C     THE ARRAY  W  IS USED AS SCRATCH SPACE.
C     ----------------------------------------------------------
C
C     +++++ TERMINATION CODES AND INTERMEDIATE PRINTING +++++
C     ----------------------------------------------------------
C     MXS  SETS A LIMIT ON THE NUMBER OF MINIMIZATION STEPS TO BE
C     TAKEN.
C
C     UPON TERMINATION  IFL  WILL BE SET ACCORDING TO
C     THE FOLLOWING CODE ...
C
C             IFL = 1 .... SUCCESSFUL TERMINATION.
C
C             IFL = 2 .... UNSUCCESSFUL TERMINATION.
C                          CONSTRAINTS CANNOT BE SATISFIED.
C                          PROBLEM IS INFEASIBLE.
C
C             IFL = 3 .... LIMIT IMPOSED BY  MXS  REACHED
C                          WITHOUT FINDING A SOLUTION.
C
C             IFL = 4 .... PROGRAM ABORTED.
C                          NUMERICAL DIFFICULTIES
C                          DUE TO ILL-CONDITIONING.
C
C             IFL = 5 .... NEQNS, NVARS, NEQC AND/OR
C                          NIQC  HAVE IMPROPER VALUES
C                          OR  E  CONTAINS A ZERO COLUMN.
C
C     IN ALL CASES THE OUTPUT PARAMETERS  X,EL1N AND RES
C     WILL CONTAIN THE VALUES WHICH THEY REACHED AT TERMINATION.
C
C     INTERMEDIATE PRINTING WILL BE TURNED OFF IF  PSW = .FALSE. .
C     ON THE OTHER HAND,  DETAILS OF EACH MINIMIZATION CYCLE
C     WILL BE PRINTED IF  PSW  IS SET TO  .TRUE.
C     ----------------------------------------------------------
C
C     +++++ REMARKS AND USER CAUTIONS +++++
C     ----------------------------------------------------------
C       1.  BEYOND SOME PRECAUTIONARY STEPS TAKEN IN
C           CERTAIN DIVISIONS, NO SPECIAL
C           OVERFLOW/UNDERFLOW PROTECTION IS PROVIDED.
C       2.  ALL TOLERANCES FOR CHECKING ZEROS AND LINEAR
C           DEPENDENCIES ARE DETERMINED FROM THE QUANTITY
C           EPS  WHICH APPEARS IN DATA DECLARATIONS IN  CL1
C           AND SEVERAL OF ITS SUBROUTINES.   EPS  CAN BE SET TO THE
C           LEAST POSITIVE NUMBER SATISFYING  (1.0 + EPS) .GT. 1.0
C           IN THE PRECISION OF ARITHMETIC BEING USED.  WITH THIS
C           SETTING,  CL1  USES AN EXTREMELY STRICT ZERO TOLERANCE.
C           FOR A MORE FORGIVING VERSION OF  CL1,  EPS MAY BE REMOVED
C           FROM THE DATA DEFINITIONS AND INCLUDED AS A USER-SPECIFIED
C           ZERO-TESTING PARAMETER IN THE ARGUMENT LIST.  IN SUCH AN
C           EVENT, IF THE PROBLEM DATA IS GIVEN TO  NDIG  SIGNIFICANT
C           DIGITS, THEN  10.0**(-NDIG)  IS A REASONABLE CHOICE
C           FOR THE VALUE OF  EPS.
C           OVERFLOW CHECKING PRIOR TO DIVISION IS
C           DONE USING THE QUANTITY  BIG,  ALSO SPECIFIED AS DATA.
C           BIG  SHOULD BE THE LARGEST REPRESENTABLE FLOATING
C           POINT NUMBER.
C       3.  THIS IS A SINGLE PRECISION VERSION OF  CL1.
C           TO CHANGE THIS CODE INTO DOUBLE PRECISION ...
C           A.  CHANGE ALL OCCURRENCES OF  - REAL -
C               DECLARATIONS TO  - DOUBLE PRECISION -
C           B.  CHANGE ALL OCCURRENCES OF  - SYSTEM ROUTINES -
C               (AS LISTED IN THE HEADING OF EACH SUBROUTINE)
C               TO THEIR CORRESPONDING DOUBLE PRECISION VERSIONS
C           C.  CHANGE ALL OCCURRENCES OF THE STRINGS  E+  AND  E-
C               TO  D+  AND  D-  RESPECTIVELY
C           D.  CHANGE ALL BASIC LINEAR ALGEBRA ROUTINES  (BLAS)
C               TO THEIR DOUBLE-PRECISION EQUIVALENTS
C           E.  BOTH  EPS  AND  BIG  WILL HAVE TO BE CHANGED
C               TO THEIR DOUBLE PRECISION EQUIVALENTS
C           F.  THE REFERENCES TO  - IFIX -  IN SUBROUTINES
C               - RESID -  AND  - GETV -  MUST BE CHANGED
C               FROM THE FORM
C                    IFIX(FLOAT(K)*UNIF(...))
C               TO THE FORM
C                    IFIX(FLOAT(K)*SNGL(UNIF(...)))
C           G.  THE REFERENCE TO  - FLOAT -  IN SUBROUTINE
C               - RESID -  MUST BE CHANGED FROM THE FORM
C                    FLOAT(...)
C               TO THE FORM
C                    DBLE(FLOAT(...))
C           H.  REMOVE THESE COMMENT CARDS  (3., 3.A.-3.H.).
C     ----------------------------------------------------------
C
C     +++++ PORTABILITY NOTE +++++
C     ----------------------------------------------------------
C     THE INTRINSIC FORTRAN FUNCTIONS WHICH HAVE BEEN USED
C     ARE NOT DECLARED IN TYPE STATEMENTS IN  CL1  OR ANY
C     OF ITS SUBROUTINES.  A VERY FEW FORTRAN COMPILERS DO
C     REQUIRE SUCH TYPING, AND THEY WILL FAIL QUITE EXPLICITLY
C     BECAUSE OF OUR OMISSION.  ON THE OTHER HAND, SOME FEW
C     OTHER FORTRAN COMPILERS DO NOT EXPECT SUCH TYPING, AND
C     THEIR MODE OF FAILURE IS MORE SUBTILE.  THESE COMPILERS
C     MAY ATTACH THE NAMES OF THE INTRINSIC FUNCTIONS TO THE
C     LIST OF SUBROUTINES TO BE OBTAINED FROM EXTERNAL SOURCES,
C     CAUSING A FAULT IN LIBRARY SEARCHING DURING PRE-EXECUTION
C     LOADING.   SINCE THERE IS NO WAY TO WIN THE GAME, WE
C     HAVE OPTED TO LOSE IT, IN THOSE FEW CASES WHERE WE MUST,
C     IN THE MOST IMMEDIATE AND OBVIOUS FASHION.
C     NOTE -- INTRINSIC FUNCTIONS ARE LISTED AMONG THE
C                  - SYSTEM ROUTINES -
C              AT THE START OF EACH SUBPROGRAM.
C
C     THE QUANTITIES  EPS  AND  BIG  WHICH ARE TO BE FOUND
C     IN DATA STATEMENTS ARE MACHINE DEPENDENT.  SEE  2.  UNDER
C     + REMARKS AND USER CAUTIONS +.
C
C     MANY OF THE COMPUTATIONS IN THIS PROGRAM ARE
C     DEFINED IN TERMS OF THE  BASIC LINEAR ALGEBRA
C     SUBROUTINES  (BLAS).  REFER TO THE APPENDIX
C     OF THIS CODE FOR FURTHER INFORMATION.
C     **********************************************************
C
      INTEGER DDX,GRDX,ICYC,IADDC,IDELC
      INTEGER PX,PTEX,RRX,TOPX,ZZX
      REAL ALPHA,AMAG,CGMAG,PEN,PENPAR
C
C     ////////////////  BEGIN PROGRAM  /////////////////////////
C
      CALL SETUP
     *                    (NEQNS,NEQC,NIQC,NVARS,DDX,GRDX,PX,PTEX,RRX,
     *                     TOPX,ZZX,ICYC,IFL,E,NER,AMAG,CGMAG,PENPAR)
C
C
   10   CONTINUE
C
          CALL NEWPEN
     *                    (IADDC,IDELC,NACT,NEQNS,NEQC,
     *                     NIQC,NVARS,IFL,E,NER,X,F,RES,
     *                     W(PTEX),ALPHA,PENPAR,INDX)
C
   20     CONTINUE
C
            CALL UPDATE
     *                    (IADDC,IDELC,NACT,NEQNS,NEQC,NIQC,NVARS,
     *                     ICYC,IFL,MXS,E,NER,X,F,RES,
     *                     W(GRDX),EL1N,PEN,PENPAR,
     *                     INDX,W(ZZX),NVARS,W(DDX),W(RRX),
     *                     W(TOPX))
            CALL MONIT
     *                    (NACT,NEQC,NIQC,NVARS,ICYC,PSW,X,
     *                     ALPHA,EL1N,PEN,PENPAR,INDX)
            CALL FINDP
     *                    (IDELC,NACT,NEQNS,NEQC,NIQC,NVARS,IFL,E,NER,
     *                     X,F,RES,W(GRDX),W(PX),EL1N,AMAG,CGMAG,
     *                     PENPAR,INDX,W(ZZX),NVARS,W(DDX),
     *                     W(RRX),W(TOPX))
            CALL STEP
     *                    (IADDC,NACT,NEQNS,NEQC,NIQC,NVARS,IFL,E,NER,
     *                     X,RES,W(GRDX),W(PX),W(PTEX),ALPHA,PENPAR,
     *                     INDX,W(TOPX))
C
            IF (IFL .EQ. 0)  GO TO 20
C
          IF (IFL .EQ. 2
     *           .AND. (CGMAG + PENPAR*AMAG) .NE. CGMAG)  GO TO 10
C
      RETURN
      END
C
C     ---------------
C     SECOND LEVEL SUBROUTINES --
C          SETUP,NEWPEN,UPDATE,MONIT,FINDP,STEP,REFINE
C     ---------------
C
      SUBROUTINE SETUP
     *                    (NEQNS,NEQC,NIQC,NVARS,DDX,GRDX,PX,PTEX,RRX,
     *                     TOPX,ZZX,ICYC,IFL,E,NER,AMAG,CGMAG,PENPAR)
C
      INTEGER DDX,GRDX,ICYC,IFL,NEQC,NEQNS,NER
      INTEGER NIQC,NVARS,PX,PTEX,RRX,TOPX,ZZX
      REAL AMAG,CGMAG,E(NER,1),PENPAR
C
C     ***************
C     CL1  VERSION.
C
C     SET UP THE PROGRAM
C     PARAMETERS AND INDICES.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C     +++++++++++++++
C
      INTEGER I,J,NCOLS,NQNP1
      REAL OCT,TMP,ZERO
C
      DATA OCT /8.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
C     ***************
C     CHECK VALIDITY OF PROBLEM DIMENSIONS
C     ***************
C
      NCOLS = NEQNS + NEQC + NIQC
      IF (NVARS .GE. 1
     *       .AND. NEQC .GE. 0
     *        .AND. NIQC .GE. 0
     *         .AND. NEQNS .GE. 0
     *          .AND. NCOLS .GE. NVARS
     *           .AND. NER .GE. NVARS)  GO TO 10
        IFL = 5
        GO TO 100
   10 CONTINUE
C
C     ***************
C     SET UP INDICES FOR THE TEMPORARY STORAGE VECTOR  W.
C     ***************
C
      NQNP1 = NEQNS + 1
      GRDX = 1
      PX = GRDX + NVARS
      PTEX = PX + NVARS
      DDX = PTEX + NCOLS
      RRX = DDX + NVARS
      ZZX = RRX + (((NVARS + 1)*(NVARS + 2))/2)
      TOPX = ZZX + NVARS*NVARS
C
C     ***************
C     AMAG  IS A ROUGH ESTIMATE OF THE NORM OF  A.
C     CGMAG  IS A ROUGH ESTIMATE OF THE NORM OF  (G,C).
C     TOGETHER THEY ARE USED TO DETERMINE WHEN THE
C     PENALTY PARAMETER IS TOO SMALL AND WHEN THE
C     RESTRICTED GRADIENT IS ZERO.
C     ***************
C
      AMAG = ZERO
      CGMAG = ZERO
      IF (1 .GT. NEQNS)  GO TO 50
        DO 40 J=1,NEQNS
          TMP = ZERO
          DO 20 I=1,NVARS
            TMP = TMP + ABS(E(I,J))
  20     CONTINUE
          IF (TMP .GT. ZERO)  GO TO 30
            IFL = 5
            GO TO 100
   30     CONTINUE
          IF (TMP .GT. AMAG)  AMAG = TMP
   40   CONTINUE
   50 CONTINUE
      IF (NQNP1 .GT. NCOLS)  GO TO 90
        DO 80 J=NQNP1,NCOLS
          TMP = ZERO
          DO 60 I=1,NVARS
            TMP = TMP + ABS(E(I,J))
  60     CONTINUE
          IF (TMP .GT. ZERO)  GO TO 70
            IFL = 5
            GO TO 100
   70     CONTINUE
          IF (TMP .GT. CGMAG)  CGMAG = TMP
   80   CONTINUE
   90 CONTINUE
C
C     ***************
C     INITIALIZE  IFL,ICYC,PENPAR
C     ***************
C
      IFL = 2
      ICYC = -1
      PENPAR = OCT
  100 CONTINUE
      RETURN
      END
      SUBROUTINE NEWPEN
     *                    (IADDC,IDELC,NACT,NEQNS,NEQC,NIQC,
     *                     NVARS,IFL,E,NER,X,F,RES,PTE,
     *                     ALPHA,PENPAR,INDX)
C
      INTEGER IADDC,IDELC,IFL,INDX(1),NACT
      INTEGER NEQC,NEQNS,NER,NIQC,NVARS
      REAL ALPHA,E(NER,1),F(1),PENPAR,PTE(1),RES(1),X(1)
C
C     ***************
C     CL1  VERSION.
C
C     BEGIN A ROUND OF MINIMIZATION STEPS
C     WITH A NEW PENALTY PARAMETER VALUE.
C     ***************
C
C     +++++++++++++++
C     BLAS  SDOT
C     +++++++++++++++
C
      INTEGER I,NCOLS
      REAL OCT,ONE,ZERO
C
      REAL SDOT
C
      DATA ZERO /0.0E+00/
      DATA ONE /1.0E+00/
      DATA OCT /8.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
C     ***************
C     SET PENALTY PARAMETER VALUE.
C     ERASE RECORD OF ACTIVE EQUATION/CONSTRAINTS.
C     ***************
C
      IF (IFL .NE. 2)  GO TO 20
        NCOLS = NEQNS + NEQC + NIQC
        IFL = 0
        NACT = 0
        IADDC = 0
        IDELC = 0
        ALPHA = ZERO
        PENPAR = PENPAR/OCT
C
C     ***************
C     INITIALIZE  INDX,RES,PTE,INDX
C     ***************
C
        DO 10 I=1,NCOLS
          RES(I) = SDOT(NVARS,E(1,I),1,X,1) - F(I)
          PTE(I) = ZERO
          INDX(I) = I
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE UPDATE
     *                     (IADDC,IDELC,NACT,NEQNS,NEQC,NIQC,NVARS,
     *                     ICYC,IFL,MXS,E,NER,X,F,RES,GRD,
     *                     EL1N,PEN,PENPAR,INDX,ZZ,
     *                     NZZR,DD,RR,W)
C
      INTEGER IADDC,IDELC,ICYC,IFL,INDX(1),MXS
      INTEGER NACT,NEQC,NEQNS,NER,NIQC,NVARS,NZZR
      REAL DD(1),E(NER,1),EL1N,F(1),GRD(1),PEN,PENPAR
      REAL RES(1),RR(1),W(1),X(1),ZZ(NZZR,1)
C
C     ***************
C     CL1  VERSION.
C
C     PREPARATION FOR NEXT MINIMIZATION STEP.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C     +++++++++++++++
C
      INTEGER NALLQ,NCOLS
      REAL ONE,ZERO
C
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
C     ***************
C     DETERMINE THE ACTIVE EQUATIONS AND ACTIVE
C     CONSTRAINTS.  COMPUTE RESIDUALS AND FUNCTION VALUE.
C     UPDATE THE  Z*D*R  DECOMPOSITION.
C     ***************
C
      NALLQ = NEQNS + NEQC
      NCOLS = NALLQ + NIQC
      IF (IFL .NE. 0)  GO TO 20
        ICYC = ICYC + 1
        IF (ICYC .LE. MXS)  GO TO 10
          IFL = 3
          GO TO 20
   10   CONTINUE
          CALL DELCOL(IADDC,IDELC,NACT,NVARS,ZZ,NZZR,DD,RR,INDX)
          CALL RESID(IADDC,NACT,NCOLS,NVARS,E,NER,X,F,RES,INDX)
          CALL ADDCOL(IADDC,IDELC,NACT,NVARS,ZZ,NZZR,DD,RR,
     *                E,NER,INDX,W)
          CALL OBJECT(IADDC,NACT,NEQNS,NALLQ,NCOLS,NVARS,E,NER,
     *                RES,GRD,EL1N,PEN,PENPAR,INDX)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE MONIT
     *                    (NACT,NEQC,NIQC,NVARS,ICYC,PSW,X,
     *                     ALPHA,EL1N,PEN,PENPAR,INDX)
C
      INTEGER ICYC,INDX(1),NACT,NEQC,NIQC,NVARS
      LOGICAL PSW
      REAL ALPHA,EL1N,PEN,PENPAR,X(1)
C
C     ***************
C     CL1  VERSION.
C
C     PRINT OUT INFORMATION AT
C     EACH STEP IF  PSW = .TRUE.
C
C     TO CHANGE THE OUTPUT MEDIUM,
C     CHANGE THE  DATA  DECLARATION
C     APPEARING BELOW.
C     ***************
C
      INTEGER I,NXOUTP
C
      DATA NXOUTP /6/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (.NOT. PSW)  GO TO 20
        WRITE (NXOUTP,6000) ICYC,ALPHA
        WRITE (NXOUTP,6010) (X(I),I=1,NVARS)
        WRITE (NXOUTP,6020)  NACT
        IF (NACT .LE. 0)  GO TO 10
          WRITE (NXOUTP,6030)
          WRITE (NXOUTP,6040) (INDX(I),I=1,NACT)
   10   CONTINUE
        WRITE (NXOUTP,6050) EL1N
        IF (NEQC + NIQC .GT. 0)  WRITE (NXOUTP,6060) PENPAR,PEN
   20 CONTINUE
      RETURN
C
 6000 FORMAT(20H0***** CYCLE NUMBER ,I5/
     *       6X,13HSTEP TAKEN = ,1PE15.7/
     *       6X,11HX-VECTOR...)
 6010 FORMAT(5X,3(1PE15.7))
 6020 FORMAT(6X,38HNUMBER OF ACTIVE EQUATIONS/CONSTRAINTS,I5)
 6030 FORMAT(6X,31HACTIVE EQUATIONS/CONSTRAINTS...)
 6040 FORMAT(5X,7I5)
 6050 FORMAT(6X,19HRESIDUAL L1 NORM = ,1PE15.7)
 6060 FORMAT(6X,20HPENALTY PARAMETER = ,1PE10.2/
     *       6X,25HPENALTY FUNCTION VALUE = ,1PE15.7)
      END
      SUBROUTINE FINDP
     *                    (IDELC,NACT,NEQNS,NEQC,NIQC,NVARS,IFL,
     *                     E,NER,X,F,RES,GRD,P,EL1N,AMAG,CGMAG,
     *                     PENPAR,INDX,ZZ,NZZR,DD,RR,W)
C
      INTEGER IDELC,IFL,INDX(1),NACT,NEQC,NEQNS,NER,NIQC,NVARS,NZZR
      REAL AMAG,CGMAG,DD(1),E(NER,1),EL1N,F(1),GRD(1),P(1),PENPAR
      REAL RES(1),RR(1),W(1),X(1),ZZ(NZZR,1)
C
C     ***************
C     CL1  VERSION.
C
C     DETERMINE DESCENT DIRECTION  P
C     (OR DISCOVER OPTIMALITY)
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C
C     BLAS  SASUM,SCOPY,SSCAL
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES   (1.0 + EPS) .GT. 1.0   IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     (ALTERNATIVELY, FOR LESS STRICT ZERO CHECKING,
C      EPS  CAN BE SET TO A USER-SPECIFIED TOLERANCE.)
C     +++++++++++++++
C
      INTEGER COEFX,I,IX,NALLQ,NALQP1,NCOLS,NQNP1,TOPX
      LOGICAL FAIL
      REAL EPS,GRDNRM,ONE,PNRM,PROD,TEST,ZERO
C
      REAL SASUM
C
C     +++++++++++++++
C     THE FOLLOWING DECLARATIONS ARE NECESSARY
C     FOR PORTABILITY WHEN  SCOPY  IS USED, AS
C     IT IS BELOW, TO FILL ARRAYS WITH A SINGLE VALUE
C     (ONE=UNITY  AND  ZERO=ZIP  IN THIS CASE).
C     +++++++++++++++
C
      REAL UNITY(1),ZIP(1)
      EQUIVALENCE (ONE,UNITY(1)),(ZERO,ZIP(1))
C
      DATA EPS /1.19209290E-07/
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IDELC = 0
      IF (IFL .NE. 0)  GO TO 90
        NALLQ = NEQNS + NEQC
        NALQP1 = NALLQ + 1
        NCOLS = NALLQ + NIQC
        NQNP1 = NEQNS + 1
        COEFX = 1
        TOPX = COEFX + NVARS
C
C     ***************
C     PROJECT THE NEGATIVE OF THE RESTRICTED GRADIENT
C     ONTO THE ORTHOGONAL COMPLEMENT OF THE SPACE
C     SPANNED BY THE ACTIVE COLUMNS.
C     ***************
C
        CALL ZDRPOC(NVARS,NACT,ZZ,NZZR,DD,GRD,P,FAIL)
        IF (.NOT.FAIL)  GO TO 5
          IFL = 4
          GO TO 90
    5   CONTINUE
        CALL SSCAL(NVARS,-ONE,P,1)
        PNRM = SASUM(NVARS,P,1)
        GRDNRM = SASUM(NVARS,GRD,1)
C
C     ***************
C     IF THE PROJECTION IS NOT ZERO,
C     IT WILL SERVE AS A DESCENT DIRECTION.
C
C     OTHERWISE FIND THE REPRESENTATION OF
C     THE RESTRICTED GRADIENT AS A LINEAR
C     COMBINATION OF THE ACTIVE COLUMNS.
C     THE COEFFICIENTS OF THE LINEAR COMBINATION
C     ARE TO BE STORED IN THE ARRAY  COEF
C     (THAT IS, IN  W(COEFX),...,W(COEFX+NACT-1)).
C     ***************
C
        IF (PNRM .GT. EPS*(AMAG*PENPAR + CGMAG))  GO TO 90
        IF (NACT .EQ. 0)  GO TO 50
          CALL ZDRGNV(NVARS,NACT,ZZ,NZZR,RR,GRD,W(COEFX),FAIL)
          IF (.NOT. FAIL)  GO TO 10
            IFL = 4
            GO TO 90
   10     CONTINUE
C
C     ***************
C     CONVERT THE COEFFICIENTS OF THE LINEAR
C     COMBINATION INTO A DESCENT DIRECTION  P ,
C     OR DETERMINE OPTIMALITY.
C
C     IF THE OPTIMALITY TEST IS NOT SATISFIED,
C     GETV  WILL INDICATE AN EQUATION/CONSTRAINT
C     TO BE DELETED FROM ACTIVITY BY THE VALUE
C     OF  IDELC.  FOR OPTIMALITY,  IDELC=0.
C     ***************
C
          CALL GETV(IDELC,NACT,NVARS,NEQNS,NALLQ,E,NER,GRD,
     *              W(COEFX),PENPAR,INDX)
          PNRM = ZERO
          IF (IDELC .EQ. 0)  GO TO 20
            CALL ZDRGIT(NVARS,NACT,ZZ,NZZR,RR,W(COEFX),P,FAIL,W(TOPX))
            IF (.NOT. FAIL)  PNRM = SASUM(NVARS,P,1)
            IF (.NOT. FAIL)  GO TO 20
              IFL = 4
              GO TO 90
   20   CONTINUE
C
C     ***************
C     IF A DESCENT DIRECTION  P  COULD HAVE BEEN FOUND,
C     IT HAS BEEN OBTAINED BY THIS POINT IN THE PROGRAM.
C
C     CHECK FOR OPTIMALITY.
C
C     PNRM  HAS BEEN SET EXACTLY ZERO
C     AFTER THE CALL TO SUBROUTINE  GETV
C     IF THE OPTIMALITY CONDITIONS ARE SATISFIED.
C     THE CHECK BELOW HAS BEEN MADE SOMEWHAT
C     COMPLICATED TO ALLOW FOR THE RARE EVENT THAT
C     THE RESTRICTED GRADIENT IS ZERO AND NO
C     COLUMNS ARE ACTIVE,  OR THAT THE  L1  NORM OF
C               (A-TRANSPOSE) * X - F
C     IS COMPUTATIONALLY ZERO.
C     (THE CALL TO THE SUBROUTINE  REFINE
C      MAY BE OMITTED, IF DESIRED.)
C     ***************
C
        IF (PNRM .LE. EPS*(AMAG*PENPAR + CGMAG))  GO TO 50
        DO 40 I=1,NEQNS
          TEST = ABS(F(I))
          DO 30 IX=1,NVARS
            PROD = ABS(E(IX,I)*X(IX))
            IF (PROD .GT. TEST)  TEST = PROD
   30     CONTINUE
          IF (ABS(RES(I)) .GT. EPS*TEST)  GO TO 90
   40   CONTINUE
   50   CONTINUE
          IFL = 1
          CALL REFINE(NACT,NEQNS,NCOLS,NVARS,IFL,E,NER,
     *                X,F,EL1N,RES,INDX,ZZ,NZZR,RR,W)
          IF (IFL .NE. 1)  GO TO 90
C
C     ***************
C     IF THE PROBLEM HAS CONSTRAINTS,
C     CHECK FEASIBILITY.
C     ***************
C
          IF (NQNP1 .GT. NCOLS)  GO TO 90
            DO 80 I=NQNP1,NCOLS
              TEST = ABS(F(I))
              DO 60 IX=1,NVARS
                PROD = ABS(E(IX,I)*X(IX))
                IF (PROD .GT. TEST)  TEST = PROD
   60         CONTINUE
              TEST = EPS*TEST
              IF (I .GT. NALLQ)  GO TO 70
                IF (ABS(RES(I)) .LE. TEST)  GO TO 80
                  IFL = 2
                  GO TO 90
   70         CONTINUE
                IF (RES(I) .GE. (-TEST))  GO TO 80
                  IFL = 2
                  GO TO 90
   80       CONTINUE
   90 CONTINUE
      RETURN
      END
      SUBROUTINE STEP
     *                   (IADDC,NACT,NEQNS,NEQC,NIQC,NVARS,IFL,
     *                    E,NER,X,RES,GRD,P,PTE,ALPHA,PENPAR,INDX,ALF)
C
      INTEGER IADDC,IFL,INDX(1),NACT,NEQC,NEQNS,NER,NIQC,NVARS
      REAL ALPHA,ALF(1),E(NER,1),GRD(1),P(1)
      REAL PENPAR,PTE(1),RES(1),X(1)
C
C     ***************
C     CL1  VERSION.
C
C     PIECEWISE LINEAR LINE SEARCH.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES ABS,SIGN
C
C     BLAS  SASUM,SAXPY,SDOT
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES   (1.0 + EPS) .GT. 1.0   IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     (ALTERNATIVELY, FOR LESS STRICT ZERO CHECKING,
C      EPS  CAN BE SET TO A USER-SPECIFIED TOLERANCE.)
C
C     BIG  IS THE LARGEST POSITIVE NUMBER
C     WHICH CAN BE REPRESENTED IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     +++++++++++++++
C
      INTEGER I,IIN,IX,JX,NACTP1,NALLQ,NCOLS,NUM,NUMNAC
      REAL BIG,DEN,EPS,GRDNRM,ONE,PNRM,PTG
      REAL RATIO,RESID,TMP,TWO,ZERO
C
      REAL SASUM,SDOT
C
      DATA BIG /1.0E+30/
      DATA EPS /1.19209290E-07/
      DATA ONE /1.0E+00/
      DATA TWO /2.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
C     ***************
C     THIS ROUTINE DETERMINES ALL OF THE RATIOS  ALF
C     OF THE FORM
C        -RES(I)/((E(.,I)-TRANSP)*P),
C              FOR  I = K+1,...,MPL
C     WHICH ARE NONNEGATIVE AND HENCE INDICATE DISTANCES
C     FROM THE POINT  X  TO BREAKPOINTS WHICH WILL
C     BE ENCOUNTERED IN TRAVEL ALONG DIRECTION  P.
C     THE INDEX VECTOR  INDX  IS REARRANGED SO THAT
C     ITS  K+1  THROUGH  NUM  COMPONENTS CORRESPOND TO
C     THESE NONNEGATIVE RATIOS.
C     THE RESULTS ARE HEAPED SO THAT THE  ALF  VALUES CAN
C     BE INSPECTED IN ORDER FROM SMALLEST TO LARGEST.
C     THE BREAKPOINT  ALPHA  GIVING THE MINIMUM OBJECTIVE
C     FUNCTION VALUE IS FOUND, AND  X  IS
C     ADJUSTED TO  X + ALPHA*P .
C
C     THE INNER PRODUCTS  (E(.,I)-TRANSPOSE)*P  ARE SAVED
C     FOR LATER USE IN UPDATING THE RESIDUAL VALUES.
C     ***************
C
      ALPHA = ZERO
      IF (IFL .NE. 0)  GO TO 90
        NALLQ = NEQNS + NEQC
        NCOLS = NALLQ + NIQC
        NACTP1 = NACT + 1
        NUM = 0
        IF (1 .GT. NACT)  GO TO 20
          DO 10 I=1,NACT
            IX = INDX(I)
            PTE(IX) = SDOT(NVARS,E(1,IX),1,P,1)
   10     CONTINUE
   20   CONTINUE
        IF (NACTP1 .LE. NCOLS)  GO TO 30
          IFL = 1
          GO TO 90
   30   CONTINUE
          DO 50 I=NACTP1,NCOLS
            IX = INDX(I)
            RESID = RES(IX)
            DEN = SDOT(NVARS,E(1,IX),1,P,1)
            PTE(IX) = DEN
            IF (SIGN(ONE,RESID) .EQ. SIGN(ONE,DEN)
     *                    .AND. RESID .NE. ZERO)  GO TO 50
              RESID = ABS(RESID)
              DEN = ABS(DEN)
              IF (DEN .GE. ONE)  GO TO 40
                IF (RESID .GE. DEN*BIG)  GO TO 50
   40         CONTINUE
                RATIO = RESID/DEN
                NUM = NUM + 1
                NUMNAC = NUM + NACT
                JX = INDX(NUMNAC)
                INDX(NUMNAC) = IX
                INDX(I) = JX
                ALF(NUM) = RATIO
   50     CONTINUE
          IF (NUM .GT. 0)  GO TO 60
            IFL = 2
            GO TO 90
   60     CONTINUE
C
C     ***************
C     HEAP THE POSITIVE RATIOS
C     ***************
C
            CALL DKHEAP(.TRUE.,NUM,INDX(NACTP1),ALF)
C
C     ***************
C     TRAVEL ALONG  P  UNTIL NO FURTHER DECREASE IN THE
C     PENALTY FUNCTION IS POSSIBLE
C     ***************
C
            IIN = NUM
            PTG = SDOT(NVARS,GRD,1,P,1)
            PNRM = SASUM(NVARS,P,1)
            GRDNRM = SASUM(NVARS,GRD,1)
            DO 70 I=1,NUM
              IX = INDX(NACTP1)
              TMP = -SIGN(ONE,RES(IX))
              IF (IX .LE. NALLQ)  TMP = TMP*TWO
              IF (IX .LE. NEQNS)  TMP = TMP*PENPAR
              PTG = PTG + TMP*PTE(IX)
              IF (PTG .GE. (-EPS*GRDNRM*PNRM))  GO TO 80
                CALL DKHEAP(.FALSE.,IIN,INDX(NACTP1),ALF)
   70       CONTINUE
            IFL = 2
            GO TO 90
   80       CONTINUE
            IADDC = NACTP1
C
C     ***************
C     ADJUST  X  TO  X + ALPHA*P
C     ***************
C
            ALPHA = ALF(1)
            CALL SAXPY(NVARS,ALPHA,P,1,X,1)
   90 CONTINUE
      RETURN
      END
      SUBROUTINE REFINE(NACT,NEQNS,NCOLS,NVARS,IFL,E,NER,
     *                  X,F,EL1N,RES,INDX,ZZ,NZZR,RR,W)
C
      INTEGER IFL,INDX(1),NACT,NCOLS,NEQNS,NER,NVARS,NZZR
      REAL E(NER,1),EL1N,F(1),RES(1),RR(1),W(1),X(1),ZZ(NZZR,1)
C
C     ***************
C     A ROUTINE FOR REFINING THE SOLUTION
C     PRODUCED BY  CL1.
C
C     (THIS ROUTINE MAY BE OMITTED IF DESIRED.)
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C
C     BLAS  SDOT
C     +++++++++++++++
C
      INTEGER I,IX
      LOGICAL FAIL
      REAL TMP,ZERO
C
      REAL SDOT
C
      DATA ZERO /0.0E+00/
C
C     /////////////// BEGIN PROGRAM ///////////////
C
      IF (NACT .EQ. 0)  GO TO 40
        DO 10 I=1,NACT
          IX = INDX(I)
          RES(I) = F(IX)
   10   CONTINUE
        CALL ZDRGIT(NVARS,NACT,ZZ,NZZR,RR,RES,X,FAIL,W)
        IF (.NOT. FAIL)  GO TO 20
          IFL = 4
          GO TO 40
   20   CONTINUE
        EL1N = ZERO
        DO 30 I=1,NCOLS
          TMP = SDOT(NVARS,E(1,I),1,X,1) - F(I)
          RES(I) = TMP
          IF (I .LE. NEQNS)  EL1N = EL1N + ABS(TMP)
   30   CONTINUE
   40 CONTINUE
      RETURN
      END
C
C     ---------------
C     THIRD LEVEL SUBROUTINES --
C          DELCOL,RESID,ADDCOL,OBJECT,GETV
C     ---------------
C
      SUBROUTINE DELCOL(IADDC,IDELC,NACT,NROW,ZZ,NZZR,DD,RR,INDX)
C
      INTEGER INDX(1),NZZR,NACT,IDELC,IADDC,NROW
      REAL DD(1),RR(1),ZZ(NZZR,1)
C
C     ***************
C     CL1  VERSION OF  DELCOL.
C
C     THIS ROUTINE ADMINISTERS THE DELETION OF THE COLUMN
C     INDICATED BY THE VALUE OF IDELC
C     FROM AN  NROW BY NACT   Z*D*R   DECOMPOSITION.
C     NOTE THAT THE VALUE OF IDELC
C     IS THE NUMBER OF A COLUMN IN THE DECOMPOSITION
C     RATHER THAN A NUMBER WHICH REFERS TO
C     A COLUMN IN THE MATRIX  E.
C     (THE  E-COLUMN  NUMBERS CORRESPONDING TO
C     THE COLUMNS OF THE FACTORIZATION ARE TO BE
C     FOUND IN   INDX(1),...,INDX(NACT) .
C     THE CONTENTS OF   INDX(NACT+1),...,INDX(IADDC)
C     INDICATE COLUMNS OF  E  WHICH ARE SLATED FOR
C     ADDITION TO THE DECOMPOSITION.)
C     THE VECTOR  INDX   IS REARRANGED BY
C     PERMUTING THE ELEMENT WHICH CORRESPONDS TO
C     THE DELETION OUT TO THE   IADDC-TH  POSITION.
C     NACT  AND  IADDC  ARE DECREASED ACCORDINGLY.
C     ***************
C
      INTEGER I,IDLP1,IXDLC
      LOGICAL FAIL
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (IDELC .EQ. 0)  GO TO 20
        IDLP1 = IDELC + 1
        IXDLC = INDX(IDELC)
        DO 10 I=IDLP1,IADDC
          INDX(I - 1) = INDX(I)
   10   CONTINUE
        INDX(IADDC) = IXDLC
        IADDC = IADDC - 1
        CALL ZDRCOU(NROW,NACT,ZZ,NZZR,DD,RR,IDELC,FAIL)
        IDELC = IXDLC
   20 CONTINUE
      RETURN
      END
      SUBROUTINE RESID(IADDC,NACT,NCOLS,NVARS,E,NER,X,F,RES,INDX)
C
      INTEGER IADDC,INDX(1),NACT,NCOLS,NER,NVARS
      REAL E(NER,1),F(1),RES(1),X(1)
C
C
C     ***************
C     COMPUTE THE RESIDUALS
C          (E(.,IX)-TRANSP)*X - F(IX)  .
C     THE RESIDUALS ARE STORED IN THE ARRAY  RES.
C     INDX  IS REARRANGED SO THAT THE ZERO RESIDUALS
C     CORRESPOND TO  INDX(1),...,INDX(IADDC)  .
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS,IFIX,FLOAT,SQRT
C
C     BLAS  SDOT
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES   (1.0 + EPS) .GT. 1.0   IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     (ALTERNATIVELY, FOR LESS STRICT ZERO CHECKING,
C      EPS  CAN BE SET TO A USER-SPECIFIED TOLERANCE.)
C     +++++++++++++++
C
      INTEGER I,IADP1,IDUMMY,IRAND,IX,J,NACTP1
      REAL EPS,PROD,TEMP,TEST,TOL,ZERO
C
      REAL SDOT,UNIF01
C
      DATA EPS /1.19209290E-07/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      TOL = EPS*SQRT(FLOAT(NVARS))
      NACTP1 = NACT + 1
      IF (1 .GT. IADDC)  GO TO 20
C
C     ***************
C     ZERO OUT ALL RESIDUALS KNOWN TO BE ZERO.
C     ***************
C
        DO 10 I=1,IADDC
          IX = INDX(I)
          RES(IX) = ZERO
   10   CONTINUE
   20 CONTINUE
C
C     ***************
C     COMPUTE THE REMAINING RESIDUALS.
C     DETECT ANY MORE RESIDUALS WHICH
C     ARE COMPUTATIONALLY ZERO, AND
C     SET THEM EXACTLY ZERO.  THEIR
C     ASSOCIATED INDICES ARE PERMUTED
C     SO THAT THEY ARE STORED IN
C     INDX(NACT+1),...,NACT(IADDC).
C
C     (A FAIRLY TIGHT ZERO CHECK IS USED.
C     IT IS FAR LESS EXPENSIVE IN RUNNING
C     TIME TO NEGLECT AN EXTRA ZERO
C     RESIDUAL THAN TO ACCEPT IT AND RISK
C     INVOKING THE ANTI-CYCLING
C     MECHANISMS IN THE PROGRAM.
C     THE ACCURACY OF THE SOLUTION AS
C     FINALLY DETERMINED IS NOT AFFECTED.)
C     ***************
C
      IADP1 = IADDC + 1
      IF (IADP1 .GT. NCOLS)  GO TO 50
        DO 40 I=IADP1,NCOLS
          IX = INDX(I)
          TEMP = SDOT(NVARS,E(1,IX),1,X,1) - F(IX)
          TEST = ABS(F(IX))
          DO 25 J=1,NVARS
            PROD = ABS(E(J,IX)*X(J))
            IF (PROD .GT. TEST)  TEST = PROD
   25     CONTINUE
          TEST = TOL*TEST
          IF (ABS(TEMP) .LE. TEST)  GO TO 30
            RES(IX) = TEMP
            GO TO 40
   30     CONTINUE
            IADDC = IADDC + 1
            INDX(I) = INDX(IADDC)
            INDX(IADDC) = IX
            RES(IX) = ZERO
   40   CONTINUE
   50 CONTINUE
C
C     ***************
C     IF ANY NEW ZERO RESIDUALS HAVE
C     BEEN FOUND, RANDOMIZE THEIR
C     ORDERING AS AN ANTI-CYCLING
C     DEVICE FOR  ADDCOL.
C     ***************
C
      IF (IADDC .LE. NACTP1)  GO TO 70
        DO 60 I=NACTP1,IADDC
          IRAND = I + IFIX(FLOAT(IADDC - I + 1)*UNIF01(0,IDUMMY))
          IX = INDX(IRAND)
          INDX(IRAND) = INDX(I)
          INDX(I) = IX
   60   CONTINUE
   70 CONTINUE
      RETURN
      END
      SUBROUTINE ADDCOL(IADDC,IDELC,NACT,NVARS,ZZ,NZZR,DD,RR,
     *                  E,NER,INDX,W)
C
      INTEGER IADDC,IDELC,INDX(1),NACT,NER,NVARS,NZZR
      REAL DD(1),E(NER,1),RR(1)
      REAL W(1),ZZ(NZZR,1)
C
C     ***************
C     CL1 VERSION OF ADDCOL.
C
C     THIS ROUTINE ADMINISTERS THE ADJUSTMENT OF THE
C     Z*D*R   DECOMPOSITION FOR ANY NEW ZERO RESIDUALS.
C     THE DATA CORRESPONDING TO THE ZERO RESIDUALS IS INDEXED
C     IN  INDX(NACT+1),...,INDX(IADDC).
C     ***************
C
C     +++++++++++++++
C     BLAS  SASUM
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES   (1.0 + EPS) .GT. 1.0   IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     (ALTERNATIVELY, FOR LESS STRICT ZERO CHECKING,
C      EPS  CAN BE SET TO A USER-SPECIFIED TOLERANCE.)
C     +++++++++++++++
C
      INTEGER I,ISTRT,IX,NACTP1,TOPX
      LOGICAL FAIL
      REAL COLNRM,EPS,PRJNRM
C
      REAL SASUM
C
      DATA EPS /1.19209290E-07/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      TOPX = NVARS + 1
      ISTRT = NACT + 1
      IF (ISTRT .GT. IADDC)  GO TO 20
C
C     ***************
C     CANDIDATES FOR ADDITION TO THE  Z*D*R
C     FACTORIZATION ARE INSPECTED IN RANDOM
C     ORDER TO HINDER CYCLING.
C     THE RANDOMIZATION WAS CARRIED OUT BY  RESID.
C
C     IF A CANDIDATE HAS JUST BEEN RELEASED
C     FROM THE FACTORIZATION OR IS DEPENDENT UPON THE
C     COLUMNS IN THE FACTORIZATION,
C     THEN IT IS OMITTED FROM ADDITION.
C
C     UPON EXIT, INDICES OF SUCH OMITTED
C     COLUMNS ARE TO BE FOUND IN
C          INDX(NACT+1),...,INDX(IADDC) .
C     ***************
C
        DO 10 I=ISTRT,IADDC
          NACTP1 = NACT + 1
          IX = INDX(I)
          CALL ZDRPOC(NVARS,NACT,ZZ,NZZR,DD,E(1,IX),W,FAIL)
          COLNRM = SASUM(NVARS,E(1,IX),1)
          PRJNRM = SASUM(NVARS,W,1)
          IF (PRJNRM .LE. EPS*COLNRM .OR. IX .EQ. IDELC)  GO TO 10
            INDX(I) = INDX(NACTP1)
            INDX(NACTP1) = IX
            CALL ZDRCIN(NVARS,NACT,ZZ,NZZR,DD,RR,E(1,IX),FAIL,W)
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE OBJECT(IADDC,NACT,NEQNS,NALLQ,NCOLS,NVARS,
     *                  E,NER,RES,GRD,EL1N,PEN,PENPAR,INDX)
C
      INTEGER IADDC,INDX(1),NACT,NALLQ,NEQNS,NCOLS,NER,NVARS
      REAL E(NER,1),EL1N,GRD(1),PEN,PENPAR,RES(1)
C
C     ***************
C     CL1 VERSION OF OBJECT.
C
C     THIS ROUTINE ADMINISTERS THE EVALUATION OF THE
C     PENALTY (OBJECTIVE) FUNCTION GIVEN THE EQUATION
C     AND CONSTRAINT RESIDUALS.  IT ALSO COMPUTES THE
C     RESTRICTED GRADIENT OF THE FUNCTION.
C
C     COLUMNS WHICH ARE NOT IN THE  Z*D*R FACTORIZATION
C     BUT WHICH ARE ASSOCIATED WITH ZERO RESIDUALS MUST
C     BE INCLUDED IN THE RESTRICTED GRADIENT WITH RANDOM
C     SIGNS AS AN ANTI-CYCLING DEVICE.
C     THE INDICES OF THESE COLUMNS ARE TO BE
C     FOUND IN  INDX(NACT+1),...,INDX(IADDC)
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS,SIGN
C
C     BLAS  SAXPY,SCOPY
C     +++++++++++++++
C
      INTEGER I,IDUMMY,IX,NACTP1
      REAL HALF,ONE,SGN,TMP,ZERO
C
      REAL UNIF01
C
C     +++++++++++++++
C     THE FOLLOWING DECLARATIONS ARE NECESSARY
C     FOR PORTABILITY WHEN  SCOPY  IS USED, AS
C     IT IS BELOW, TO FILL ARRAYS WITH A SINGLE
C     VALUE  (ZERO=ZIP  IN THIS CASE).
C     +++++++++++++++
C
      REAL ZIP(1)
      EQUIVALENCE (ZERO,ZIP(1))
C
      DATA HALF /0.5E+00/
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      NACTP1 = NACT + 1
      EL1N = ZERO
      PEN = ZERO
      CALL SCOPY(NVARS,ZIP,0,GRD,1)
      IF (NACTP1 .GT. NCOLS)  GO TO 30
        DO 20 I=NACTP1,NCOLS
          IX = INDX(I)
          TMP = RES(IX)
          SGN = SIGN(ONE,TMP)
          IF (I .LE. IADDC)  SGN = SIGN(ONE,HALF - UNIF01(0,IDUMMY))
          IF (IX .GT. NALLQ .AND. SGN .GT. ZERO)  GO TO 20
            TMP = ABS(TMP)
            IF (IX .GT. NEQNS)  GO TO 10
              EL1N = EL1N + TMP
              TMP = TMP*PENPAR
              SGN = SGN*PENPAR
   10       CONTINUE
            PEN = PEN + TMP
            CALL SAXPY(NVARS,SGN,E(1,IX),1,GRD,1)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE GETV(IDELC,NACT,NVARS,NEQNS,NALLQ,E,NER,
     *                GRD,COEF,PENPAR,INDX)
C
      INTEGER IDELC,INDX(1),NACT,NALLQ,NEQNS,NER,NVARS
      REAL COEF(1),E(NER,1),GRD(1),PENPAR
C
C     ***************
C     CL1  VERSION.
C
C     SET UP THE RIGHT-HAND-SIDE VECTOR
C     (AND STORE IN THE ARRAY  COEF)
C     FOR THE LINEAR PROBLEM WHICH DETERMINES
C     A DESCENT DIRECTION  P  IN THE CASE WHERE
C     THE PROJECTION OF THE RESTRICTED GRADIENT IS ZERO.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS,FLOAT,IFIX,SIGN
C
C     BLAS  SAXPY
C
C     EPS  IS THE SMALLEST POSITIVE NUMBER WHICH
C     SATISFIES   (1.0 + EPS) .GT. 1.0   IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     (ALTERNATIVELY, FOR LESS STRICT ZERO CHECKING,
C      EPS  CAN BE SET TO A USER-SPECIFIED TOLERANCE.)
C     +++++++++++++++
C
      INTEGER I,IDUMMY,IRAND,IX
      REAL CF,EPS,ONE,OPE,S,TMP,TMPSAV,ZERO
C
      REAL UNIF01
C
      DATA EPS /1.19209290E-07/
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
C     ***************
C     FIND THE MOST OUT-OF-KILTER
C     COEFFICIENT.  BEGIN INSPECTING
C     THE COEFFICIENTS AT A RANDOM INDEX
C     TO HINDER CYCLING.  SET  COEF
C     TO ZERO ON THE FLY.
C     ***************
C
      OPE = ONE + EPS
      IDELC = 0
      TMPSAV = ZERO
      IF (1 .GT. NACT)  GO TO 40
        IRAND = IFIX(FLOAT(NACT)*UNIF01(0,IDUMMY))
        DO 30 I=1,NACT
          IRAND = IRAND + 1
          IF (IRAND .GT. NACT)  IRAND = 1
          IX = INDX(IRAND)
          CF = COEF(IRAND)
          COEF(IRAND) = ZERO
          IF (IX .GT. NALLQ)  GO TO 10
            IF (IX .LE. NEQNS)  CF = CF/PENPAR
            TMP = OPE - ABS(CF)
            GO TO 20
   10     CONTINUE
            TMP = CF + EPS
   20    CONTINUE
         IF (TMP .GE. TMPSAV)  GO TO 30
           IDELC = IRAND
           S = SIGN(ONE,CF)
           TMPSAV = TMP
   30   CONTINUE
C
C     ***************
C     IF NO COEFFICIENTS ARE OUT OF KILTER,
C     THEN RETURN.  OTHERWISE SET A
C     VALUE IN AN APPROPRIATE COMPONENT
C     (INDICATED BY  IDELC)  OF  COEF
C     AND ADJUST THE RESTRICTED GRADIENT
C     IF NECESSARY.
C     ***************
C
        IF (IDELC .EQ. 0)  GO TO 40
          COEF(IDELC) = ONE
          IX = INDX(IDELC)
          IF (IX .GT. NALLQ)  GO TO 40
            COEF(IDELC) = -S
            TMP = -S
            IF (IX .LE. NEQNS)  TMP = TMP*PENPAR
            CALL SAXPY(NVARS,TMP,E(1,IX),1,GRD,1)
   40 CONTINUE
      RETURN
      END
C
C     ---------------
C     FOURTH LEVEL SUBROUTINES --
C               DKHEAP,UNIF01,ZDRCIN,ZDRCOU,
C               ZDRGIT,ZDRGNV,ZDRPOC
C     ---------------
C
      SUBROUTINE DKHEAP(MAKE,IR,INDX,ARAY)
C
      INTEGER INDX(1),IR
      LOGICAL MAKE
      REAL ARAY(1)
C
C     ***************
C     AN ADAPTATION OF D. E. KNUTH,S HEAPING
C     ROUTINES (SEE VOLUME 3 OF
C          THE ART OF COMPUTER PROGRAMMING  ).
C     IF  MAKE  IS  .TRUE.,  THE FULL HEAP BUILDING
C     PROCESS IS CARRIED OUT ON
C          ARAY(1),...,ARAY(IR) ,
C     AND THE VALUE OF  IR  IS UNCHANGED.
C     IF  MAKE  IS  .FALSE.,  ONE STEP OF THE SORTING
C     PROCESS IS CARRIED OUT TO PROVIDE THE NEXT
C     ELEMENT OF  ARAY  IN ORDER,  AND THE VARIABLE
C     IR  IS DECREASED BY ONE.  THE INTERRUPTION OF THE
C     SORTING PHASE IS BUILT IN VIA THE FLAG  ONCE.
C     INDX  IS AN INDEX VECTOR ASSOCIATED WITH
C     ARAY  WHICH MUST BE REARRANGED IN PARALLEL
C     WITH IT.
C     ***************
C
      INTEGER I,IL,IT,J
      LOGICAL ONCE
      REAL T
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (IR .GT. 1)  GO TO 5
        IF (.NOT.MAKE)  IR = 0
        RETURN
    5 CONTINUE
C
C     ***************
C     TEST WHETHER OR NOT THE INITIAL
C     HEAP IS TO BE BUILT
C     ***************
C
      IL = 1
      IF (MAKE)  IL = (IR/2) + 1
      ONCE = .FALSE.
C
C     ***************
C     THE LOOP BEGINS HERE
C     ***************
C
   10 CONTINUE
      IF (IL .GT. 1)  GO TO 20
C
C     ***************
C     THE SORTING PHASE USES THIS BRANCH
C     ***************
C
        IF (MAKE .OR. ONCE)  RETURN
        ONCE = .TRUE.
        IT = INDX(IR)
        T = ARAY(IR)
        INDX(IR) = INDX(1)
        ARAY(IR) = ARAY(1)
        IR = IR - 1
        IF (IR .GT. 1)  GO TO 30
          INDX(1) = IT
          ARAY(1) = T
          RETURN
   20 CONTINUE
C
C     ***************
C     THE HEAP-BUILDING PHASE USES THIS BRANCH
C     ***************
C
        IL = IL - 1
        IT = INDX(IL)
        T = ARAY(IL)
   30 CONTINUE
C
C     ***************
C     THE REMAINING STATEMENTS ARE COMMON
C     TO BOTH PHASES AND EMBODY THE
C     HEAP-RECTIFYING (SIFTING) SECTION
C     ***************
C
      J = IL
   40 CONTINUE
      I = J
      J = 2*J
      IF (J - IR)  50, 60, 70
   50 CONTINUE
      IF (ARAY(J) .LE. ARAY(J + 1))  GO TO 60
        J = J + 1
   60 CONTINUE
      IF (T .LE. ARAY(J))  GO TO 70
        INDX(I) = INDX(J)
        ARAY(I) = ARAY(J)
        GO TO 40
   70 CONTINUE
      INDX(I) = IT
      ARAY(I) = T
      GO TO 10
      END
      REAL FUNCTION UNIF01(ISEED,IX)
C
      INTEGER ISEED,IX,IX0
C
      DATA IX0/2/
C
C     +++++++++++++++
C     SYSTEM ROUTINES  FLOAT,MOD
C     +++++++++++++++
C
C     --------------------------------------------------------------
C     --------------------------------------------------------------
C
C     *****PURPOSE-
C     THIS FUNCTION RETURNS A PSEUDO-RANDOM NUMBER DISTRIBUTED
C     UNIFORMLY IN THE INTERVAL (0,1).
C
C     *****PARAMETER DESCRIPTION-
C     ON INPUT-
C
C     ISEED,  IF IT IS NONZERO MODULO 9973, BECOMES THE
C          NEW SEED, I.E. IT REPLACES THE INTERNALLY STORED
C          VALUE OF IX0.  ON MACHINES WHERE FORTRAN VARIABLES
C          RETAIN THEIR VALUES BETWEEN CALLS, THE INTERNALLY
C          STORED VALUE IF IX0 IS THE VALUE ASSIGNED TO  IX  IN
C          THE PREVIOUS INVOCATION OF  UNIF01.  OTHERWISE -- AND
C          IN THE FIRST CALL TO  UNIF01 --  IX0=2.
C
C     ON OUTPUT-
C
C     IX IS THE NEXT INTEGER IN A PSEUDO-RANDOM SEQUENCE OF
C          INTEGERS BETWEEN  1  AND  9972  AND IS GENERATED FROM ITS
C          PREDECESSOR  IX0  (I.E.  FROM  ISEED,  IF  ISEED  IS NONZERO
C          MODULO 9973).  IX  IS THE VALUE WHICH  ISEED  SHOULD HAVE
C          IN THE NEXT INVOCATION OF  UNIF01  TO GET THE NEXT
C          PSEUDO-RANDOM NUMBER.  THE CALLER WILL OFTEN PASS THE
C          SAME VARIABLE FOR  ISEED  AS FOR  IX,
C          E.G.  X = UNIF01(IX,IX).
C
C     *****APPLICATION AND USAGE RESTRICTIONS-
C     UNIF01  SHOULD ONLY BE USED WHEN PORTABILITY IS IMPORTANT AND A
C     COURSE RANDOM NUMBER GENERATOR SUFFICES.  APPLICATIONS REQUIRING
C     A FINE, HIGH PRECISON GENERATOR SHOULD USE ONE WITH A MUCH
C     LARGER MODULUS.
C
C     *****ALGORITHM NOTES-
C     UNIF01 WILL RUN ON ANY MACHINE HAVING AT LEAST 20 BITS OF AC-
C     CURACY FOR FIXED-POINT ARITHMITIC.  IT IS BASED ON A GENERATOR
C     RECOMMENDED IN (3), WHICH PASSES THE SPECTRAL TEST WITH FLYING
C     COLORS -- SEE (1) AND (2).
C
C     REFERENCES-
C     (1) HOAGLIN, D.C. (1976), THEORETICAL PROPERTIES OF CONGRUENTIAL
C     RANDOM-NUMBER GENERATORS-  AN EMPIRICAL VIEW,
C     MEMORANDUM NS-340, DEPT. OF STATISTICS, HARVARD UNIV.
C
C     (2) KNUTH, D.E. (1969), THE ART OF COMPUTER PROGRAMMING, VOL. 2
C     (SEMINUMERICAL ALGORITHMS), ADDISON-WESLEY, READING, MASS.
C
C     (3) SMITH, C.S. (1971), MULTIPLICATIVE PSEUDO-RANDOM NUMBER
C     GENERATORS WITH PRIME MODULUS, J. ASSOC. COMPUT. MACH. 18,
C     PP. 586-593.
C
C     *****GENERAL-
C
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH
C     SUPPORTED BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS
C     MCS-7600324, DCR75-10143, 76-14311DSS, AND MCS76-11989.
C
C     PERMISSION FOR THE USE OF  UNIF01  IN  CL1  WAS
C     GENEROUSLY GIVEN BY  V. KLEMA  AND  D. HOAGLIN.
C
C     --------------------------------------------------------------
C     --------------------------------------------------------------
C
      IF (ISEED .EQ. 0) GO TO 10
      IX = MOD(ISEED, 9973)
      IF (IX .NE. 0) IX0 = IX
C  ***
C  IN ORDER THAT ALL FIXED-POINT CALCULATIONS REQUIRE ONLY 20 BIT
C  ARITHMETIC, WE USE TWO CALLS TO  MOD  TO COMPUTE
C  IX0 = MOD(3432*IX0, 9973).
C  ***
 10   IX0 = MOD(52*MOD(66*IX0, 9973), 9973)
      IX = IX0
      UNIF01 = FLOAT(IX0)/9973.0E+00
      RETURN
      END
      SUBROUTINE ZDRCIN(N,K,ZZ,NZZR,DD,RR,COL,FAIL,W)
C
      INTEGER K,N,NZZR
      LOGICAL FAIL
      REAL COL(1),DD(1),RR(1),W(1),ZZ(NZZR,1)
C
C     ***************
C     PREPARED BY RICHARD BARTELS
C     THE UNIVERSITY OF WATERLOO
C     COMPUTER SCIENCE DEPARTMENT
C     LATEST UPDATE .... 30 NOVEMBER, 1979.
C
C     GIVEN THE FACTORIZATION
C
C          ZZ*DD*RR
C
C     OF SOME  N BY K  MATRIX
C
C       (0 .LE. K .LT. N)
C         (N .GE. 1),
C
C     WHERE
C
C          (ZZ-TRANSP)*(ZZ) = (DD-INV),
C          DD  IS DIAGONAL AND NONSINGULAR,
C     AND
C          RR  HAS ZEROS BELOW THE DIAGONAL,
C
C     AND GIVEN A  (K+1)TH  COLUMN
C     TO BE ADDEDTO THE ORIGINAL MATRIX,
C     THIS PROGRAM UPDATES  ZZ,DD AND RR.
C
C     THE VALUE OF  K  IS INCREASED BY ONE.
C
C     W  IS A SCRATCH ARRAY.
C
C     USE IS MADE OF ROUTINES FROM THE LIBRARY
C     OF BASIC LINEAR ALGEBRA SUBROUTINES (BLAS).
C
C     PARAMETERS...
C
C                     INPUT/
C       NAME   TYPE   OUTPUT/   SUB-    DESCRIPTION
C                     SCRATCH  SCRIPTS
C       -------------------------------------------
C       N      INT.      I              NUMBER OF ROWS
C
C       K      INT.     I/O             NUMBER OF COLUMNS
C
C       ZZ     REAL     I/O       2     SCALED ORTHOGONAL
C                                       MATRIX
C
C       NZZR   INT.      I              ROW DIMENSION OF ZZ
C
C       DD     REAL     I/O       1     DIAGONAL SCALING
C                                       MATRIX (DIAGONAL
C                                       ELEMENTS ONLY)
C
C       RR     REAL     I/O       1     RIGHT-TRIANGULAR
C                                       MATRIX IN COMPACT FORM.
C
C       COL    REAL      I        1     COLUMN TO BE
C                                       ADDED TO  RR
C
C       FAIL   LOG.      O             .TRUE.  IF  K,N
C                                       ARE IMPROPER
C
C       W      REAL     SCR       1     WORKSPACE
C       -------------------------------------------
C
C     THE  I-TH  SEGMENT OF THE ARRAY  RR  IS  N-I+2 SPACES
C     LONG AND CONTAINS  1  WORK SPACE FOLLOWED BY THE
C     K-I+1  ELEMENTS OF ROW  I  FOLLOWED BY  N-K
C     SCRATCH SPACES.
C     ***************
C
C     +++++++++++++++
C     BLAS  SCOPY,SDOT,SROTM,SROTMG
C     +++++++++++++++
C
      INTEGER I,J,JDEL,KP1,KP2
      REAL DI,ONE,PARAM(5),WI,ZERO
C
      REAL SDOT
C
C     +++++++++++++++
C     THE FOLLOWING DECLARATIONS ARE NECESSARY
C     FOR PORTABILITY WHEN  SCOPY  IS USED, AS
C     IT IS BELOW, TO FILL ARRAYS WITH A SINGLE VALUE
C     (ONE=UNITY  AND  ZERO=ZIP  IN THIS CASE).
C     +++++++++++++++
C
      REAL UNITY(1),ZIP(1)
      EQUIVALENCE (ONE,UNITY(1)),(ZERO,ZIP(1))
C
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (K .GE. 0 .AND. K .LT. N .AND. N .LE. NZZR)  GO TO 5
        FAIL = .TRUE.
        RETURN
    5 CONTINUE
      IF (K .GT. 0)  GO TO 20
C
C     ***************
C     FOR THE SPECIAL CASE THAT THE
C     FACTORIZATION WAS VACUOUS,
C     RESET THE ARRAYS  ZZ AND DD
C     TO REPRESENT THE  IDENTITY.
C     ***************
C
        K = 0
        CALL SCOPY(N,UNITY,0,DD,1)
        CALL SCOPY(N*N,ZIP,0,ZZ,1)
        DO 10 I=1,N
          ZZ(I,I) = ONE
   10   CONTINUE
   20 CONTINUE
      KP1 = K + 1
      KP2 = K + 2
C
C     ***************
C     TRANSFORM THE INCOMING COLUMN,
C     AND STORE THE RESULT IN  W.
C     ***************
C
      DO 30 I=1,N
        W(I) = SDOT(N,ZZ(1,I),1,COL,1)
   30 CONTINUE
C
C     ***************
C     ZERO OUT THE SPIKE WHICH WOULD RESULT FROM
C     STORING  W  IN  RR.   UPDATE  ZZ  AND  DD.
C     ***************
C
      IF (KP2 .GT. N)  GO TO 50
        DO 40 I=KP2,N
          DI = DD(I)
          WI = W(I)
          CALL SROTMG(DD(KP1),DI,W(KP1),WI,PARAM)
          W(I) = WI
          DD(I) = DI
          CALL SROTM(N,ZZ(1,KP1),1,ZZ(1,I),1,PARAM)
   40   CONTINUE
   50 CONTINUE
C
C     ***************
C     STORE THE NEW COLUMN, WHICH IS STILL
C     IN THE ARRAY  W,  INTO  RR.
C     ***************
C
      J = KP2
      JDEL = N
      DO 60 I=1,KP1
        RR(J) = W(I)
        J = J + JDEL
        JDEL = JDEL - 1
   60 CONTINUE
      K = KP1
      FAIL = .FALSE.
      RETURN
      END
      SUBROUTINE ZDRCOU(N,K,ZZ,NZZR,DD,RR,IC,FAIL)
C
      INTEGER IC,K,N,NZZR
      LOGICAL FAIL
      REAL DD(1),RR(1),ZZ(NZZR,1)
C
C     ***************
C     PREPARED BY RICHARD BARTELS
C     THE UNIVERSITY OF WATERLOO
C     COMPUTER SCIENCE DEPARTMENT
C     LATEST UPDATE .... 30 NOVEMBER, 1979.
C
C     GIVEN THE FACTORIZATION
C
C          ZZ*DD*RR
C
C     OF SOME  N BY K  MATRIX
C
C       (1 .LE. K .LE. N)
C          (N .GE. 1),
C
C     WHERE
C
C          (ZZ-TRANSP)*(ZZ) = (DD-INV),
C          DD  IS DIAGONAL AND NONSINGULAR,
C     AND
C          RR  HAS ZEROS BELOW THE DIAGONAL,
C
C     AND GIVEN THE INDEX  IC  OF A COLUMN
C     TO BE REMOVED  (1 .LE. IC .LE. K),
C     THIS PROGRAM UPDATES  ZZ,DD AND RR .
C
C     THE VALUE OF  K  IS DECREASED BY ONE, AND
C     THE COLUMN ORDERING IN  RR  IS CHANGED.
C
C     USE IS MADE OF ROUTINES FROM THE LIBRARY
C     OF BASIC LINEAR ALGEBRA SUBROUTINES (BLAS).
C
C     PARAMETERS...
C
C       NAME   TYPE    INPUT/   SUB-    DESCRIPTION
C                      OUTPUT  SCRIPTS
C       -------------------------------------------
C       N      INT.      I              NUMBER OF ROWS
C
C       K      INT.     I/O             NUMBER OF COLUMNS
C
C       ZZ     REAL     I/O       2     SCALED ORTHOGONAL
C                                       MATRIX
C
C       NZZR   INT.      I              ROW DIMENSION OF ZZ
C
C       DD     REAL     I/O       1     DIAGONAL SCALING
C                                       MATRIX (DIAGONAL
C                                       ELEMENTS ONLY)
C
C       RR     REAL     I/O       1     RIGHT-TRIANGULAR
C                                       MATRIX IN COMPACT FORM.
C
C       IC     INT.      I              INDEX OF COLUMN
C                                       TO BE REMOVED
C
C       FAIL   LOG.      O              .TRUE.  IF  K,N,IC
C                                       ARE IMPROPER
C       -------------------------------------------
C
C     THE  I-TH  SEGMENT OF THE ARRAY  RR  IS  N-I+2 SPACES
C     LONG AND CONTAINS  1  WORK SPACE FOLLOWED BY THE
C     K-I+1  ELEMENTS OF ROW  I  FOLLOWED BY  N-K
C     SCRATCH SPACES.
C     ***************
C
C     +++++++++++++++
C     BLAS  SROTM,SROTMG
C     +++++++++++++++
C
      INTEGER I,IM1,J,JEND,JINC,JSTRT,KM1,LSTRT
      REAL DI,PARAM(5),RJ
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (K .GE. 1 .AND. K .LE. N .AND. N .LE. NZZR)  GO TO 5
        FAIL = .TRUE.
        RETURN
    5 CONTINUE
      KM1 = K - 1
C
C     ***************
C     SPECIAL CASES ARE HANDLED FIRST.
C     1.  K=1 AND THE FACTORIZATION BECOMES NULL.
C     2.  IC=K AND THE UPDATING IS TRIVIAL.
C     ***************
C
      IF (K .GT. 1)  GO TO 10
        K = 0
        FAIL = .FALSE.
        RETURN
   10 CONTINUE
      IF (IC .LT. K)  GO TO 20
        K = KM1
        FAIL = .FALSE.
        RETURN
   20 CONTINUE
C
C     ***************
C     GENERAL UPDATING STEP.
C     THE COLUMN TO BE DELETED MUST BE PERMUTED
C     TO THE RIGHT, AND SUBDIAGONAL ELEMENTS
C     WHICH RESULT IN  RR  HAVE TO BE
C     TRANSFORMED TO ZERO.
C     ***************
C
      JSTRT = IC + 1
      JEND = K
      JINC = N
      DO 50 I=1,K
C
C     ***************
C     PERMUTATION OF THE  I-TH  ROW OF RR.
C     ***************
C
        DO 30 J=JSTRT,JEND
          RR(J) = RR(J + 1)
   30   CONTINUE
        IF (I .LE. IC)  GO TO 40
C
C     ***************
C     TRANSFORMATION OF THE CURRENT AND LAST
C     ROWS  (I AND I-1)  OF RR  AS WELL AS
C     CORRESPONDING CHANGES TO  ZZ AND DD.
C
C     THE EXTRA VARIABLES  DI  AND  RJ
C     ARE USED TO AVOID AN ERROR MESSAGE
C     FROM THE  PFORT VERIFIER, AND THEY
C     MAY BE REMOVED, IF DESIRED, SO THAT
C     THE CALL TO  SROTMG  WOULD BE
C
C     CALL SROTMG(DD(IM1),DD(I),RR(LSTRT),RR(JSTRT),PARAM)
C
C     ***************
C
          IM1 = I - 1
          DI = DD(I)
          RJ = RR(JSTRT)
          CALL SROTMG(DD(IM1),DI,RR(LSTRT),RJ,PARAM)
          RR(JSTRT) = RJ
          DD(I) = DI
          CALL SROTM(JEND - JSTRT + 1,RR(LSTRT + 1),1,
     *               RR(JSTRT + 1),1,PARAM)
          CALL SROTM(N,ZZ(1,IM1),1,ZZ(1,I),1,PARAM)
          JSTRT = JSTRT + 1
   40   CONTINUE
C
C     ***************
C     INDEX UPDATING
C     ***************
C
        LSTRT = JSTRT
        JSTRT = JSTRT + JINC
        JEND = JEND + JINC
        JINC = JINC - 1
   50 CONTINUE
      K = KM1
      FAIL = .FALSE.
      RETURN
      END
      SUBROUTINE ZDRGIT(N,K,ZZ,NZZR,RR,GV,SOL,FAIL,W)
C
      INTEGER K,N,NZZR
      LOGICAL FAIL
      REAL GV(1),RR(1),W(1),SOL(1),ZZ(NZZR,1)
C
C     ***************
C     PREPARED BY RICHARD BARTELS
C     THE UNIVERSITY OF WATERLOO
C     COMPUTER SCIENCE DEPARTMENT
C     LATEST UPDATE .... 30 NOVEMBER, 1979.
C
C     GIVEN THE FACTORIZATION
C
C          ZZ*DD*RR
C
C     OF SOME  N BY K  MATRIX
C
C       (1 .LE. K .LE. N)
C         (N .GE. 1),
C
C     WHERE
C
C          (ZZ-TRANSP)*(ZZ) = (DD-INV),
C          DD  IS DIAGONAL AND NONSINGULAR,
C     AND
C          RR  HAS ZEROS BELOW THE DIAGONAL,
C
C     AND GIVEN AN ARBITRARY VECTOR  GV  OF
C     APPROPRIATE DIMENSION, THIS ROUTINE FINDS THE
C     VECTOR  SOL  SATISFYING THE UNDERDETERMINED SYSTEM
C
C          (ZZ*DD*RR-TRANSP.)*(SOL) = (GV).
C
C     THAT IS,
C
C          (SOL) = ((ZZ*DD*RR)-GEN.INV.-TRANSP.)*(GV).
C
C     THE ARRAY  DD  IS NOT NEEDED BY  ZDRGIT.
C
C     USE IS MADE OF ROUTINES FROM THE LIBRARY
C     OF BASIC LINEAR ALGEBRA SUBROUTINES (BLAS).
C
C     W  IS A SCRATCH ARRAY.
C
C     PARAMETERS...
C
C                     INPUT/
C       NAME   TYPE   OUTPUT/   SUB-    DESCRIPTION
C                     SCRATCH  SCRIPTS
C       -------------------------------------------
C       N      INT.      I              NUMBER OF ROWS
C
C       K      INT.     I/O             NUMBER OF COLUMNS
C
C       ZZ     REAL     I/O       2     SCALED ORTHOGONAL
C                                       MATRIX
C
C       NZZR   INT.      I              ROW DIMENSION OF ZZ
C
C       RR     REAL     I/O       1     RIGHT-TRIANGULAR
C                                       MATRIX IN COMPACT FORM.
C
C       GV     REAL      I        1     GIVEN VECTOR
C
C       SOL    REAL      O        1     SOLUTION
C
C       FAIL   LOG.      O              .TRUE. IF  N,K
C                                       ARE IMPROPER, OR IF
C                                       RR  IS SINGULAR
C
C       W      REAL     SCR       1     WORKSPACE
C       -------------------------------------------
C
C     THE  I-TH  SEGMENT OF THE ARRAY  RR  IS  N-I+2 SPACES
C     LONG AND CONTAINS  1  WORK SPACE FOLLOWED BY THE
C     K-I+1  ELEMENTS OF ROW  I  FOLLOWED BY  N-K
C     SCRATCH SPACES.
C
C     IF  GV  AND  SOL  ARE DIMENSIONED TO THE
C     MAXIMUM OF  N  AND  K , THEN THE SAME
C     STORAGE ARRAY MAY BE USED FOR BOTH OF
C     THESE VECTORS.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C
C     BLAS  SAXPY,SCOPY
C
C     BIG  IS THE LARGEST POSITIVE NUMBER
C     WHICH CAN BE REPRESENTED IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     +++++++++++++++
C
      INTEGER I,J,JDEL
      REAL BIG,WI,ONE,RRJ,ZERO
C
C     +++++++++++++++
C     THE FOLLOWING DECLARATIONS ARE NECESSARY
C     FOR PORTABILITY WHEN  SCOPY  IS USED, AS
C     IT IS BELOW, TO FILL ARRAYS WITH A SINGLE VALUE
C     (ZERO=ZIP  IN THIS CASE).
C     +++++++++++++++
C
      REAL ZIP(1)
      EQUIVALENCE (ZERO,ZIP(1))
C
      DATA BIG /1.0E+30/
      DATA ONE /1.0E+00/
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (K .GE. 1 .AND. K .LE. N .AND. N .LE. NZZR)  GO TO 10
        FAIL = .TRUE.
        RETURN
   10 CONTINUE
C
C     ***************
C     FIRST SOLVE  (RR-TRANSP.)*(W) = (GV)
C     ***************
C
      CALL SCOPY(K,GV,1,W,1)
      J = 2
      JDEL = N + 1
      DO 30 I=1,K
        RRJ = RR(J)
        WI = W(I)
        IF (ABS(RRJ) .GE. ONE)  GO TO 20
        IF (ABS(WI) .LT. ABS(RRJ)*BIG)  GO TO 20
          FAIL = .TRUE.
          RETURN
   20   CONTINUE
        W(I) = WI/RRJ
        IF (I.LT.K) CALL SAXPY(K-I,(-W(I)),RR(J+1),1,W(I+1),1)
        J = J + JDEL
        JDEL = JDEL - 1
   30 CONTINUE
C
C     ***************
C     NOW  (SOL) = (ZZ)*(W)
C     ***************
C
      CALL SCOPY(N,ZIP,0,SOL,1)
      DO 40 I=1,K
        CALL SAXPY(N,W(I),ZZ(1,I),1,SOL,1)
   40 CONTINUE
      FAIL = .FALSE.
      RETURN
      END
      SUBROUTINE ZDRGNV(N,K,ZZ,NZZR,RR,GV,SOL,FAIL)
C
      INTEGER K,N,NZZR
      LOGICAL FAIL
      REAL GV(1),RR(1),SOL(1),ZZ(NZZR,1)
C
C     ***************
C     PREPARED BY RICHARD BARTELS
C     THE UNIVERSITY OF WATERLOO
C     COMPUTER SCIENCE DEPARTMENT
C     LATEST UPDATE .... 30 NOVEMBER, 1979.
C
C     GIVEN THE FACTORIZATION
C
C          ZZ*DD*RR
C
C     OF SOME  N BY K  MATRIX
C
C       (1 .LE. K .LE. N)
C         (N .GE. 1),
C
C     WHERE
C
C          (ZZ-TRANSP)*(ZZ) = (DD-INV),
C          DD  IS DIAGONAL AND NONSINGULAR,
C     AND
C          RR  HAS ZEROS BELOW THE DIAGONAL,
C
C     AND GIVEN AN ARBITRARY VECTOR  GV  OF
C     APPROPRIATE DIMENSION, THIS ROUTINE FINDS THE
C     VECTOR  SOL  GIVEN BY
C
C          (SOL) = ((ZZ*DD*RR)-GEN.INV.)*(GV),
C
C     WHICH REPRESENTS THE LEAST SQUARES PROBLEM
C
C        (ZZ*DD*RR)*(SOL) = (GV).
C
C     THE ARRAY  DD  IS NOT NEEDED BY  ZDRGNV.
C
C     USE IS MADE OF ROUTINES FROM THE LIBRARY
C     OF BASIC LINEAR ALGEBRA SUBROUTINES (BLAS).
C
C     PARAMETERS...
C
C       NAME   TYPE   INPUT/    SUB-    DESCRIPTION
C                     OUTPUT/  SCRIPTS
C       -------------------------------------------
C       N      INT.      I              NUMBER OF ROWS
C
C       K      INT.     I/O             NUMBER OF COLUMNS
C
C       ZZ     REAL     I/O       2     SCALED ORTHOGONAL
C                                       MATRIX
C
C       NZZR   INT.      I              ROW DIMENSION OF ZZ
C
C       RR     REAL     I/O       1     RIGHT-TRIANGULAR
C                                       MATRIX IN COMPACT FORM.
C
C       GV     REAL      I        1     GIVEN VECTOR
C
C       SOL    REAL      O        1     SOLUTION
C
C       FAIL   LOG.      O              .TRUE. IF  N,K
C                                       ARE IMPROPER, OR IF
C                                       RR  IS SINGULAR
C       -------------------------------------------
C
C     THE  I-TH  SEGMENT OF THE ARRAY  RR  IS  N-I+2 SPACES
C     LONG AND CONTAINS  1  WORK SPACE FOLLOWED BY THE
C     K-I+1  ELEMENTS OF ROW  I  FOLLOWED BY  N-K
C     SCRATCH SPACES.
C     ***************
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C
C     BLAS  SDOT
C
C     BIG  IS THE LARGEST POSITIVE NUMBER
C     WHICH CAN BE REPRESENTED IN THE
C     PRECISION OF THE ARITHMETIC BEING USED.
C     +++++++++++++++
C
      INTEGER I,IX,J,JDEL
      REAL BIG,ONE,TDEN,TNUM
C
      REAL SDOT
C
      DATA BIG /1.0E+30/
      DATA ONE /1.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      IF (K .GE. 1 .AND. K .LE. N .AND. N .LE. NZZR)  GO TO 10
        FAIL = .TRUE.
        RETURN
   10 CONTINUE
C
C     ***************
C     FORM   (V) = (ZZ(1)-TRANSP)*(GV),   WHERE  ZZ(1)
C     IS THE MATRIX OF THE FIRST  K  COLUMNS OF  ZZ
C
C     V  CAN BE STORED IN THE ARRAY  SOL.
C     ***************
C
      DO 20 I=1,K
        SOL(I) = SDOT(N,ZZ(1,I),1,GV,1)
   20 CONTINUE
C
C     ***************
C     BACKSOLVE THE SYSTEM
C          (RR)*(SOL) = (V)
C     FOR THE VECTOR  SOL
C
C     NOTE THAT  SOL  AND  V
C     ARE STORED IN THE SAME ARRAY.
C     ***************
C
      J = (((N + 1)*(N + 2) - (N - K + 3)*(N - K + 2))/2) + 2
      JDEL = N - K + 3
      DO 40 IX=1,K
        I = K - IX + 1
        TDEN = RR(J)
        TNUM = SOL(I)
        IF (IX.GT.1) TNUM = TNUM-SDOT(IX-1,RR(J+1),1,SOL(I+1),1)
        IF (ABS(TDEN) .GE. ONE)  GO TO 30
        IF (ABS(TNUM) .LT. ABS(TDEN)*BIG)  GO TO 30
          FAIL = .TRUE.
          RETURN
   30   CONTINUE
        SOL(I) = TNUM/TDEN
        J = J - JDEL
        JDEL = JDEL + 1
   40 CONTINUE
      FAIL = .FALSE.
      RETURN
      END
      SUBROUTINE ZDRPOC(N,K,ZZ,NZZR,DD,GV,POC,FAIL)
C
      INTEGER K,N,NZZR
      LOGICAL FAIL
      REAL DD(1),POC(1),GV(1),ZZ(NZZR,1)
C
C     ***************
C     PREPARED BY RICHARD BARTELS
C     THE UNIVERSITY OF WATERLOO
C     COMPUTER SCIENCE DEPARTMENT
C     LATEST UPDATE .... 30 NOVEMBER, 1979.
C
C     ZZ IS AN  N BY N  (N .GE. 1)  SCALED
C     ORTHOGONAL MATRIX.  DD  CONTAINS THE
C     DIAGONAL ELEMENTS OF A DIAGONAL SCALING
C     MATRIX.  GV  IS A GIVEN VECTOR OF LENGTH  N.
C
C     WE HAVE
C
C          (ZZ-TRANSP.)*(ZZ) = (DD-INV.)
C
C     AND
C
C               ZZ*DD*RR = MAT
C
C     FOR SOME  N BY K  (0 .LE. K .LE. N)
C     MATRIX  RR  WITH ZEROS BELOW THE DIAGONAL
C     AND SOME GIVEN MATRIX  MAT.  (NIETHER  RR
C     NOR  MAT  ARE NEEDED BY  ZDRPOC.)
C
C     THEN
C
C    (PROJ(OC)) = (ZZ(2))*(DD(2))*(ZZ(2)-TRANSP.)
C
C     IS THE (ORTHOGONAL) PROJECTOR ON THE
C     COMPLEMENT OF THE RANGE SPACE OF  MAT,
C     WHERE  ZZ(2)  REPRESENTS THE LAST  N-K
C     COLUMNS OF  ZZ  AND  DD(2)  REPRESENTS THE
C     LOWER-RIGHT-HAND  N-K  ORDER SUBMATRIX OF  DD.
C
C     ZDRPOC  PRODUCES THE VECTOR
C
C               POC = (PROJ(OC))*GV .
C
C     USE IS MADE OF ROUTINES FROM THE LIBRARY
C     OF BASIC LINEAR ALGEBRA SUBROUTINES (BLAS).
C
C     PARAMETERS...
C
C                     INPUT/
C       NAME   TYPE   OUTPUT/   SUB-    DESCRIPTION
C                     SCRATCH  SCRIPTS
C       -------------------------------------------
C       N      INT.      I              ORDER OF  ZZ,DD
C
C       K      INT.     I/O             NUMBER OF COLUMNS
C                                       OF  ZZ  DEFINING
C                                       RANGE OF  MAT
C
C       ZZ     REAL     I/O       2     SCALED ORTHOGONAL
C                                       MATRIX
C
C       NZZR   INT.      I              ROW DIMENSION OF ZZ
C
C       DD     REAL     I/O       1     DIAGONAL SCALING
C                                       MATRIX (DIAGONAL
C                                       ELEMENTS ONLY)
C
C       GV     REAL      I        1     VECTOR TO BE PROJECTED
C
C       POC    REAL      O        1     PROJECTION
C
C       FAIL   LOG.      O              .TRUE.  IF  N,K
C                                       ARE IMPROPER
C
C       -------------------------------------------
C
C     ***************
C
C     +++++++++++++++
C     BLAS  SAXPY,SCOPY,SDOT
C     +++++++++++++++
C
      INTEGER I,KP1
      REAL WI,ZERO
C
      REAL SDOT
C
C     +++++++++++++++
C     THE FOLLOWING DECLARATIONS ARE NECESSARY
C     FOR PORTABILITY WHEN  SCOPY  IS USED, AS
C     IT IS BELOW, TO FILL ARRAYS WITH A SINGLE VALUE
C     (ZERO=ZIP  IN THIS CASE).
C     +++++++++++++++
C
      REAL ZIP(1)
      EQUIVALENCE (ZERO,ZIP(1))
C
      DATA ZERO /0.0E+00/
C
C     /////////////////  BEGIN PROGRAM  //////////////////
C
      KP1 = K + 1
      IF (K .GE. 0 .AND. K .LE. N
     *       .AND. N .GE. 1 .AND. N .LE. NZZR)  GO TO 5
        FAIL = .TRUE.
        RETURN
    5 CONTINUE
      IF (K .GT. 0)  GO TO 10
C
C     ***************
C     CASE 1 ... ZZ(2)=ZZ  (K=0)
C     ***************
C
        CALL SCOPY(N,GV,1,POC,1)
        FAIL = .FALSE.
        RETURN
   10 CONTINUE
      IF (K .LT. N)  GO TO 20
C
C     ***************
C     CASE 2 ... ZZ(2) IS VACUOUS  (K=N)
C     ***************
C
        CALL SCOPY(N,ZIP,0,POC,1)
        FAIL = .FALSE.
        RETURN
   20 CONTINUE
C
C     ***************
C     CASE 3 ... ZZ(2)  IS INTERMEDIATE
C     BETWEEN THE OTHER TWO CASES
C     (0 .LT. K .LT. N)
C     ***************
C
      CALL SCOPY(N,ZIP,0,POC,1)
      DO 30 I=KP1,N
        WI = SDOT(N,ZZ(1,I),1,GV,1)*DD(I)
        CALL SAXPY(N,WI,ZZ(1,I),1,POC,1)
   30 CONTINUE
      FAIL = .FALSE.
      RETURN
      END
C
C     ---------------------------------------------
C                  ** APPENDIX **
C
C     BASIC LINEAR ALGEBRA SUBROUTINES (BLAS) --
C               SASUM,SAXPY,SCOPY,SDOT,
C               SROTM,SROTMG,SSCAL
C
C     THE  BLAS  ARE DESCRIBED IN
C
C              BASIC LINEAR ALGEBRA SUBPROGRAMS
C                    FOR FORTRAN USAGE
C
C                           BY
C              C. L. LAWSON, R. J. HANSON,
C              D. R. KINCAID, F. T. KROGH.
C              ACM TRANS. ON MATH SOFTWARE
C              VOL. 5, NO. 3, 308-325 (1979)
C
C              THE  BLAS  ARE TAKEN AS
C              PRIMITIVES.  THEIR CODE IS NOT
C              PART OF THE DEFINITION OF  CL1.
C              THEY ARE INCLUDED FOR REFERENCE
C              ONLY.  PLEASE REFER TO THE ABOVE
C              PUBLICATION FOR CURRENT INFORMATION.
C
C     IF THE  BLAS  ARE AVAILABLE INDEPENDANTLY
C     ON THE TARGET MACHINE, THIS APPENDIX SHOULD
C     BE REMOVED FROM THE CODE FOR  CL1.
C     ---------------------------------------------
C
      REAL FUNCTION SASUM(N,SX,INCX)
C
      INTEGER INCX,N
      REAL SX(1)
C
C     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
C
      INTEGER I,M,MP1,NS
      SASUM = 0.0E+00
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = N - (N/6)*6
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     *  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
      INTEGER INCX,INCY,N
      REAL SA,SX(1),SY(1)
C
C     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
C
      INTEGER I,IX,IY,M,MP1,NS
C
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
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
        SY(IY) = SY(IY) + SA*SX(IX)
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
   20 M = N - (N/4)*4
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE  SCOPY(N,SX,INCX,SY,INCY)
C
      INTEGER INCX,INCY,N
      REAL SX(1),SY(1)
C
C     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY.
C
      INTEGER I,IX,IY,M,MP1,NS
C
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
        SY(IY) = SX(IX)
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
   20 M = N - (N/7)*7
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
      REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
C
      INTEGER INCX,INCY,N
      REAL SX(1),SY(1)
C
C     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
C
      INTEGER I,IX,IY,M,MP1,NS
C
      SDOT = 0.0E+00
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = N - (N/5)*5
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     *   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
C
C        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      SUBROUTINE SROTM (N,SX,INCX,SY,INCY,SPARAM)
C
      INTEGER INCX,INCY,N
      REAL SPARAM(5),SX(1),SY(1)
C
C     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
C
C     (SX(1)     SX(N))
C     (      ...      )
C     (SY(1)     SY(N))
C
C     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
      INTEGER I,KX,KY,NSTEPS
      REAL SFLAG,SH11,SH12,SH21,SH22,TWO,W,Z,ZERO
C
      DATA ZERO,TWO /0.0E+00,2.0E+00/
C
      SFLAG=SPARAM(1)
      IF(N .LE. 0 .OR.(SFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF(SFLAG) 50,10,30
   10          CONTINUE
               SH12=SPARAM(4)
               SH21=SPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W+Z*SH12
                    SY(I)=W*SH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               SH11=SPARAM(2)
               SH22=SPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z
                    SY(I)=-W+SH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               SH11=SPARAM(2)
               SH12=SPARAM(4)
               SH21=SPARAM(3)
               SH22=SPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z*SH12
                    SY(I)=W*SH21+Z*SH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
C
          IF(SFLAG)120,80,100
   80     CONTINUE
          SH12=SPARAM(4)
          SH21=SPARAM(3)
               DO 90 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W+Z*SH12
               SY(KY)=W*SH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          SH11=SPARAM(2)
          SH22=SPARAM(5)
               DO 110 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z
               SY(KY)=-W+SH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          SH11=SPARAM(2)
          SH12=SPARAM(4)
          SH21=SPARAM(3)
          SH22=SPARAM(5)
               DO 130 I=1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z*SH12
               SY(KY)=W*SH21+Z*SH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
          END
      SUBROUTINE SROTMG (SD1,SD2,SX1,SY1,SPARAM)
C
      REAL SD1,SD2,SPARAM(5),SX1,SY1
C
C     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H
C     WHICH ZEROS THE SECOND COMPONENT OF THE 2-VECTOR
C         (SQRT(SD1)*SX1,SQRT(SD2)*SY2)**T.
C     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
C     +++++++++++++++
C     SYSTEM ROUTINES  ABS
C     +++++++++++++++
C
      INTEGER IFLAG,IGO
      REAL GAM,GAMSQ,ONE,RGAM,RGAMSQ,SFLAG,SH11,SH12
      REAL SH21,SH22,SP1,SP2,SQ1,SQ2,STEMP,SU,TWO,ZERO
C
      DATA ZERO,ONE,TWO /0.0E+00,1.0E+00,2.0E+00/,IFLAG/1/
      DATA GAM/4096.0E+00/
      DATA GAMSQ /1.678E+07/
      DATA RGAM /2.441E-04/
      DATA RGAMSQ /5.960E-08/
C
      IF(.NOT. SD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-SX1..
          GO TO 60
   10 CONTINUE
C     CASE-SD1-NONNEGATIVE
      SP2=SD2*SY1
      IF(.NOT. SP2 .EQ. ZERO) GO TO 20
          SFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      SP1=SD1*SX1
      SQ2=SP2*SY1
      SQ1=SP1*SX1
C
      IF(.NOT. ABS(SQ1) .GT. ABS(SQ2)) GO TO 40
          SH21=-SY1/SX1
          SH12=SP2/SP1
C
          SU=ONE-SH12*SH21
C
          IF(.NOT. SU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   30     CONTINUE
               SFLAG=ZERO
               SD1=SD1/SU
               SD2=SD2/SU
               SX1=SX1*SU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. SQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   50     CONTINUE
               SFLAG=ONE
               SH11=SP1/SP2
               SH22=SX1/SY1
               SU=ONE+SH11*SH22
               STEMP=SD2/SU
               SD2=SD1/SU
               SD1=STEMP
               SX1=SY1*SU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
          SFLAG=-ONE
          SH11=ZERO
          SH12=ZERO
          SH21=ZERO
          SH22=ZERO
C
          SD1=ZERO
          SD2=ZERO
          SX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. SFLAG .GE. ZERO) GO TO 90
C
          IF(.NOT. SFLAG .EQ. ZERO) GO TO 80
          SH11=ONE
          SH22=ONE
          SFLAG=-ONE
          GO TO 90
   80     CONTINUE
          SH21=-ONE
          SH12=ONE
          SFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. SD1 .LE. RGAMSQ) GO TO 130
               IF(.NOT. IFLAG.EQ.1) GO TO 105
C
C                   RECOMPUTE RESCALING PARAMETERS
C                   MORE ACCURATELY..
C
                    RGAM = ONE/GAM
                    GAMSQ = GAM**2
                    RGAMSQ = RGAM**2
                    IFLAG = 2
  105          CONTINUE
               IF(SD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               SD1=SD1*GAMSQ
               SX1=SX1*RGAM
               SH11=SH11*RGAM
               SH12=SH12*RGAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. SD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               SD1=SD1*RGAMSQ
               SX1=SX1*GAM
               SH11=SH11*GAM
               SH12=SH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. ABS(SD2) .LE. RGAMSQ) GO TO 190
               IF(SD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               SD2=SD2*GAMSQ
               SH21=SH21*RGAM
               SH22=SH22*RGAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. ABS(SD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               SD2=SD2*RGAMSQ
               SH21=SH21*GAM
               SH22=SH22*GAM
          GO TO 200
  220 CONTINUE
          IF(SFLAG)250,230,240
  230     CONTINUE
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               GO TO 260
  240     CONTINUE
               SPARAM(2)=SH11
               SPARAM(5)=SH22
               GO TO 260
  250     CONTINUE
               SPARAM(2)=SH11
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               SPARAM(5)=SH22
  260 CONTINUE
          SPARAM(1)=SFLAG
          RETURN
      END
      SUBROUTINE  SSCAL(N,SA,SX,INCX)
C
      INTEGER INCX,N
      REAL SA,SX(1)
C
C     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
C
      INTEGER I,M,MP1,NS
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = N - (N/5)*5
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END

