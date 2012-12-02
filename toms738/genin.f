      PROGRAM GENIN
C
C        *****  This is the driver for the general base programs.
C        *****  More accurately, it is a driver skeleton.
C        *****  There is a default set of test integrals
C        *****  in TESTF.  The user can replace TESTF with
C        *****  another subroutine called TESTF containing
C        *****  integrals of interest to him.
C
C   This version :  14 Aug 1992
C
C   This program tests the accuracy of numerical integration
C   using the low-discrepancy binary sequences of
C   H. Niederreiter (1988) as implemented in INLO, GOLO, and
C   related programs.  Various possible test integrals are
C   provided by the function TESTF.
C
C   Interactive input and output go through the Fortran units
C   denoted by * ;  at most installations, * corresponds to
C   the user's terminal.  For a prime-power base, Fortran unit 1
C   is used to read field arithmetic tables ;  for any base, it
C   is used to read irreducible polynomials.
C   Fortran unit 2 is used to save the output.
C
C   Our VAX implementation does not measure elapsed time.
C   It can be modified, in a system-dependent way, to do so.
C
C   GENIN and its associated subroutines are listed below.
C   An asterisk denotes a subroutine also used by GFARIT
C   (to set up the field-arithmetic tables gftabs.txt) and
C   by GFPLYS (to set up the irreducible polynomials in
C   irrtabs.txt).  
C
C       GENIN
C          INLO
C          CALCC
C          CALCV          %
C          CHARAC *       %
C          SETFLD *       %
C          PLYMUL *       %
C          GOLO
C          TESTF          %
C
C      A percent sign above denotes a routine also used
C      in the set of programs tailored for base 2.
C
C
C     Both the general-base and base-2 programs assume that
C     your computer's word length is 31 bits, excluding sign.
C     If this is not the case, modify the parameter NBITS
C     throughout the PARAMETER statements in this set of
C     programs accordingly.
C
      INTEGER MAXDIM, OUTUNT, MAXBAS, READY
      PARAMETER (MAXDIM=12, OUTUNT=2, MAXBAS = 13)
C
C   The parameter MAXDIM gives the maximum dimension that will
C   be used.  OUTUNT is the unit to save the output.
C   MAXBAS gives the maximum asymptotically-optimal base
C   used up to MAXDIM.  The user enters an appropriate value
C   of READY depending on whether or not the required files
C   indicated above have been set up.
C
      INTEGER I, NUM, DIMEN, SEQLEN, SKIP, STEP, BASE, CHARAC
      INTEGER OPTBAS(2:MAXDIM), PBASE, POWER(2:MAXBAS)
      REAL QUASI(MAXDIM), TESTF, EXACTF
      DOUBLE PRECISION SUM
C
C     The DATA statement below gives the asymptotically-optimal
C     base for the respective dimension.
C
      DATA  (OPTBAS(I), I = 2,MAXDIM) / 2,3,3,5,7,7,9,9,11,11,13 /
C
C
C     The data statement below gives values used in a possible
C     warm-up calculation.
C
      DATA  (POWER(I), I = 2,MAXBAS) / 12,8,8,6,6,6,4,4,4,4,4,4 /
C
C     There are theoretical reasons to believe that BASE ** e,
C     where e is defined for example in Bratley, Fox, and
C     Niederreiter (1991), would be a good choice.  However,
C     we don't want to come anywhere near the largest possible
C     machine-representable integer; hence, the compromise
C     exponents above.  Note: subject to this conditon,
C     it can't hurt to take an exponent greater than e, because
C     warm-up skipping of initial values is done implicitly
C     in O(1) time.  The maximum value of e for a fixed dimension
C     s grows like log s.  We allow some "fat" for the implicit
C     constant in our choice of POWER.
C
      WRITE (*,*) ' This is program GENIN'
C
      WRITE(*,*) ' If you wish to use the base 2, the '
      WRITE(*,*) ' alternative set of programs tailored for '
      WRITE(*,*) ' base 2 will run much faster. '
C
      WRITE (*,*) ' If the files gftabs.txt and irrtab.txt '
      WRITE (*,*) ' have not already been set up, then '
      WRITE (*,*) ' consult the guide to the programs for '
      WRITE (*,*) ' how to do so, set up those files per '
      WRITE (*,*) ' the guide, and [temporarily] quit this '
      WRITE (*,*) ' set of programs by entering the integer 0 '
      WRITE (*,*) ' below; otherwise, enter any other integer. '
      WRITE (*,*) ' ENTER an appropriate integer. '
      READ  (*,*)   READY
      IF (READY .EQ. 0) THEN
          WRITE (*,*) ' Set up gftabs.txt and irrtab.txt now. '
          WRITE (*,*) ' Exit GENIN. '
          STOP
      ENDIF
C
      WRITE (*,*)  'If the number of bit per word, excluding sign'
      WRITE (*,*)  'is not 31, enter the integer 0 below '
      WRITE (*,*)  'and fix the parameter NBITS everywhere; '
      WRITE (*,*)  'otherwise, enter a positive integer below.'
      WRITE (*,*)  'ENTER an appropriate integer. '
      READ  (*,*)   READY
      IF (READY .EQ. 0) THEN
          WRITE(*,*) 'Fix NBITS.'
          WRITE(*,*) 'Exit GENIN.'
          STOP
      ENDIF
C
      WRITE (*,*) ' Output file name is OUTFIL.DAT'
      OPEN (unit = OUTUNT, file = 'OUTFIL.DAT', status = 'UNKNOWN')
C
C     *****  OPEN statements are system-dependent.
C     *****  Therefore the statement above may have to be
C     *****  modified for use at your computer center.
C
C
    5 WRITE (*,*) ' Choose a test integral (1 to 4) or 0 to quit :'
      READ (*,*) NUM
      IF (NUM.LE.0) THEN
        WRITE (*,*) ' End of program GENIN'
        CLOSE (OUTUNT)
        STOP
      ENDIF
      IF (NUM.GT.4) THEN
        WRITE (*,*) ' No such test integral'
        GOTO 5
      ENDIF
C
C       *****  Each test integral is parameterized by
C       *****  its dimension.
C
   10 WRITE (*,*) ' Enter dimension :'
      READ (*,*) DIMEN
      IF (DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' Dimension may not exceed', MAXDIM
        GOTO 10
      ENDIF
C
   12 WRITE (*,*) ' Choose a prime or prime-power base .'
      WRITE (*,*) ' The asymptotically-optimal base '
      WRITE (*,*) ' for this dimension is OPTBAS(DIMEN) = '
      WRITE (*,*)   OPTBAS(DIMEN)
      WRITE (*,*) ' This base may not be empirically optimal. '
      WRITE (*,*) ' Enter BASE: '
C
      READ (*,*) BASE
      IF (CHARAC(BASE) .EQ. 0) THEN
        WRITE (*,*) ' Out of range or bad value :  try again'
        GOTO 12
      ENDIF
C
C
C        *****  The sequence length is the number of
C        *****  quasi-random points used to estimate the
C        *****  integral, excluding warm-up.
C        *****  The number of initial quasi-random points
C        *****  deleted during warm-up is given by SKIP,
C        *****  chosen below.
C
C        *****  Some users may wish to rewrite
C        *****  the driver to test a [heuristic] "convergence"
C        *****  criterion, stopping the generation of points
C        *****  when it is passed or when a specified number of
C        *****  points have been generated  -- whichever occurs
C        *****  first.
C
   15 WRITE (*,*) ' Choose sequence length :'
      WRITE (*,*) ' A power of the base is recommended; e.g.,  '
      WRITE (*,*) BASE ** POWER(BASE)
      WRITE (*,*) BASE ** ((POWER(BASE) + 1))
      WRITE (*,*) BASE ** ((POWER(BASE) + 2)) 
      WRITE (*,*) BASE ** ((POWER(BASE) + 3))
      WRITE (*,*) ' Enter SEQLEN '
      READ (*,*) SEQLEN
      IF (SEQLEN.LE.0) THEN
        WRITE (*,*) ' Length must be strictly positive'
        GOTO 15
      ENDIF
C
   20 WRITE (*,*) ' Choose the number of values to skip.'
      WRITE (*,*) ' One possibility is given by the heuristic '
      WRITE (*,*) ' formula SKIP = BASE ** POWER(BASE) '
      WRITE (*,*) ' when BASE <= MAXBAS, = 10000 otherwise '
      IF (BASE .LE. MAXBAS) THEN
          SKIP = BASE ** POWER(BASE)
        ELSE
          SKIP = 10000
      ENDIF
      WRITE (*,*) ' Numerically, this equals ', SKIP
      WRITE (8,*) ' Enter SKIP (not necessarily the value above) : '            C
      READ (*,*) SKIP
      IF (SKIP.LT.0) THEN
        WRITE (*,*) ' Number must be nonnegative'
        GOTO 20
      ENDIF
C
C
C
      CALL INLO (DIMEN, BASE, SKIP)
      WRITE (*,*) ' GENIN :  Initialization complete'
C
C Write title and the exact value of the integral
C
      WRITE (OUTUNT,27) NUM
   27 FORMAT (/,' Test integral ',I2)
      WRITE (OUTUNT,28) DIMEN, BASE, SEQLEN, SKIP
   28 FORMAT (/,' Dimension ',I6,',    Base ', I9,
     1 /,' Sequence ',I8,',    Skipped ',I7)
      WRITE (OUTUNT,29) EXACTF(NUM, DIMEN)
   29 FORMAT (/,' Correct value is ',G16.7)
      WRITE (OUTUNT,30)
   30 FORMAT(/,'      Iteration     Estimate of integral',/)
C
C Now estimate the integral
C
      WRITE (*,*)  ' The odd-looking iteration numbers '
      WRITE (*,*)  ' in the output are powers of the base.  '
C
      PBASE = BASE
      SUM = 0
      STEP = 500
      DO 100 I = 1, SEQLEN
        CALL GOLO (QUASI)
        SUM = SUM + TESTF(NUM, DIMEN, QUASI)
        IF (MOD(I,STEP).EQ.0) THEN
          WRITE (OUTUNT,99) I, SUM/I
        ENDIF
        IF (MOD(I,PBASE) .EQ. 0)  THEN
              WRITE (OUTUNT,99) I, SUM/I
              PBASE = PBASE * BASE
C
C         This finds the next power of the base.
C         There is reason to believe that convergenence
C         properties of the sequence of estimates is
C         better along the subsequence corrsponding to
C         powers of the base.
C
        ENDIF
   99     FORMAT (I12,G24.7)
          IF (I .EQ. 5000) STEP = 1000
          IF (I .EQ. 10000) STEP = 5000
  100 CONTINUE
C
      WRITE (*,*) ' GENIN :  iteration ends'
      GOTO 5
C
      END
C
C     ***** end of PROGRAM GENIN
      SUBROUTINE INLO (DIM, BASE, SKIP)
C
C   This version :  12 February 1992
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This subroutine calculates the values of Niederreiter's
C   C(I,J,R) and performs other initialization necessary
C   before calling GOLO.
C
C INPUT :
C   DIMEN - The dimension of the sequence to be generated.
C      {In the argument of INLO, DIMEN is called DIM because
C      DIMEN is subsequently passed via COMMON and is called
C      DIMEN there.}
C
C   BASE  - The prime or prime-power base to be used.
C   SKIP  - The number of values to throw away at the beginning
C           of the sequence.
C
C OUTPUT :
C   To GOLO, labelled common /COMM/.
C
C USES :
C   Calls CALCC to calculate the values of CJ.
C   Calls SETFLD to set up field arithmetic tables.
C   Calls CHARAC to check that base is a prime or prime-power
C     in the range we can handle.
C
C   -------------------------------------------------------------
C
C
C   This segment defines the common block /COMM/ and some
C   associated parameters.  These are for use in the general base
C   version of the generator.
C
      INTEGER MAXDIM, MAXFIG, NBITS
      PARAMETER (MAXDIM=12, MAXFIG=20, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   MAXFIG is the maximum number of base-Q digits we can handle.
C   MAXINT is the largest fixed point integer we can represent.
C   NBITS is the number of bits in a fixed-point integer, not
C     counting the sign.
C     ***** NBITS is machine dependent *****
C
      INTEGER C(MAXDIM, MAXFIG, 0:MAXFIG-1)
      INTEGER COUNT(0:MAXFIG-1), D(MAXDIM, MAXFIG)
      INTEGER NEXTQ(MAXDIM), QPOW(MAXFIG)
      INTEGER DIMEN, NFIGS
      REAL    RECIP
      COMMON  /COMM/ C, COUNT, D, NEXTQ, QPOW, DIMEN, NFIGS, RECIP
      SAVE    /COMM/
C
C   The common block /COMM/ :
C     C     - Contains the values of Niederreiter's C(I,J,R)
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of base-Q digits.  COUNT(R)
C             is the same as Niederreiter's AR(N) (page 54)
C             except that N is implicit.
C     D     - The values of D(I,J).
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP.
C     QPOW  - To speed things up a bit. QPOW(I) = Q ** (NFIGS-I).
C     DIMEN   - The dimension of the sequence to be generated
C     NFIGS - The number of base-Q digits we are using.
C     RECIP - 1.0 / (Q ** NFIGS)
C
C   Array C of the common block is set up by subroutine CALCC.
C   The other items in the common block are set up by INLO.
C
C   ------------------------------------------------------------
C
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication, and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, NQ, R, DIM, SKIP , BASE, CHARAC
      REAL    TEMP
C
      DIMEN = DIM
C
C     This assignment just relabels the variable
C     for subsequent use.
C
      IF (DIMEN.LE.0 .OR. DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' INLO :  Bad dimension'
        STOP
      ENDIF
C
C
      IF (CHARAC(BASE) .EQ. 0) THEN
        WRITE (*,*) ' INLO : Base not prime power or out of range'
        STOP
      ENDIF
C
      CALL SETFLD (BASE)
C
C   Calculate how many figures to use in base Q = BASE
C
      TEMP = LOG(2.0 ** NBITS - 1)/LOG(REAL(Q))
      NFIGS = MIN(MAXFIG, INT(TEMP))
C
      CALL CALCC
C
C   Set up RECIP and QPOW(I) = Q ** (NFIGS-I)
C
      RECIP = 1.0 / (Q ** NFIGS)
      QPOW(NFIGS) = 1
      DO 5 I = NFIGS-1, 1, -1
    5   QPOW(I) = Q * QPOW(I+1)
C
C   Initialize COUNT
C
      I = SKIP
      DO 10 R = 0, NFIGS-1
        COUNT(R) = MOD(I, Q)
        I = I / Q
   10 CONTINUE
      IF (I .NE. 0) THEN
        WRITE (*,*) ' INLO :  SKIP too long'
        STOP
      ENDIF
C
C   Initialize D
C
      DO 20 I = 1, DIMEN
        DO 20 J = 1, NFIGS
   20     D(I,J) = 0
C
      DO 50 R = 0, NFIGS-1
        IF (COUNT(R) .NE. 0) THEN
          DO 40 I = 1, DIMEN
            DO 30 J = 1, NFIGS
              D(I,J) = ADD (D(I,J), MUL (C(I,J,R), COUNT(R)))
   30       CONTINUE
   40     CONTINUE
        ENDIF
   50 CONTINUE
C
C   Initialize NEXTQ
C
      DO 70 I = 1, DIMEN
        NQ = 0
        DO 60 J = 1, NFIGS
          NQ = NQ + D(I,J) * QPOW(J)
   60   CONTINUE
        NEXTQ(I) = NQ
   70 CONTINUE
C
      RETURN
      END
C
C     *****  end of SUBROUTINE INLO
      SUBROUTINE CALCC
C
C    This version : 12 February 1992
C
C    See the general comments on implementing Niederreiter's
C    low-discrepancy sequences.
C
C    This routine calculates the values of the constants C(I,J,R).
C    As far as possible, we use Niederreiter's notation.
C    We calculate the values of C for each I in turn.
C    When all the values of C have been calculated, we return
C    this array to the calling program.
C
C    Irreducible polynomials are read from file 'irrtabs.txt'
C    through Fortran unit 1.  These polys must have been put on
C    the file beforehand by GFPOLYS.  Unit 1 is closed before
C    entry to CALCC and after returning from the subroutine.
C
C    Thanks to Michael Baudin for pointing out that MAXE should
C    be increased from 5 to 7, since one of the irreducible polynomials
C    that must be stored in PX has degree 7, 07 June 2010.
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication, and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
C
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C   -------------------------------------------------------------
C
C
C   This segment defines the common block /COMM/ and some
C   associated parameters.  These are for use in the general base
C   version of the generator.
C
      INTEGER MAXDIM, MAXFIG, NBITS
      PARAMETER (MAXDIM=12, MAXFIG=20, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   MAXFIG is the maximum number of base-Q digits we can handle.
C   MAXINT is the largest fixed point integer we can represent.
C   NBITS is the number of bits in a fixed-point integer, not
C     counting the sign.
C     ***** NBITS is machine dependent *****
C
      INTEGER C(MAXDIM, MAXFIG, 0:MAXFIG-1)
      INTEGER COUNT(0:MAXFIG-1), D(MAXDIM, MAXFIG)
      INTEGER NEXTQ(MAXDIM), QPOW(MAXFIG)
      INTEGER DIMEN, NFIGS
      REAL    RECIP
      COMMON  /COMM/ C, COUNT, D, NEXTQ, QPOW, DIMEN, NFIGS, RECIP
      SAVE    /COMM/
C
C   The common block /COMM/ :
C     C     - Contains the values of Niederreiter's C(I,J,R)
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of base-Q digits.  COUNT(R)
C             is the same as Niederreiter's AR(N) (page 54)
C             except that N is implicit.
C     D     - The values of D(I,J).
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP.
C     QPOW  - To speed things up a bit. QPOW(I) = Q ** (NFIGS-I).
C     DIMEN   - The dimension of the sequence to be generated
C     NFIGS - The number of base-Q digits we are using.
C     RECIP - 1.0 / (Q ** NFIGS)
C
C   Array C of the common block is set up by subroutine CALCC.
C   The other items in the common block are set up by INLO.
C
C   ------------------------------------------------------------
C
C
C
C
C
C   MAXE   - We need MAXDIM irreducible polynomials over GF(Q).
C            MAXE is the highest degree among these.
C   MAXV   - The maximum index used in V.
C   GFUNIT - The unit number to read field tables.
C   NPOLS  - The number of precalculated irreducible polys.
C
      INTEGER MAXE, MAXV, GFUNIT, NPOLS
C     PARAMETER (MAXE=5,  GFUNIT=1, NPOLS=25)
      PARAMETER (MAXE=7,  GFUNIT=1, NPOLS=25)
      PARAMETER (MAXV = MAXFIG + MAXE)
C
C INPUT :
C   DIMEN, the number of dimensions in use, and NFIGS, the number
C   of base-Q figures to be used, are passed in through common COMM.
C   Necessary field arithmetic tables are passed through common
C   FIELD.
C
C OUTPUT
C   The array C is returned to the calling program.
C
C USES
C   Subroutine CALCV is used for the auxiliary calculation
C   of the values V needed to get the Cs.
C
      INTEGER PX(-1:MAXE), B(-1:MAXDEG)
      INTEGER V(0:MAXV)
      INTEGER E, I, J, R, U
C
C   Prepare to read irreducible polynomials on unit 1.
C
      OPEN (UNIT=GFUNIT, FILE='irrtabs.txt', STATUS='old')
C
C      *****  OPEN statements are system dependent
C
   10 READ (GFUNIT, 900, END=500) I
  900 FORMAT (20I3)
      IF (I .NE. Q) THEN
        DO 20 J = 1, NPOLS
   20     READ (GFUNIT, 900)
        GOTO 10
      ENDIF
C
      DO 1000 I = 1, DIMEN
C
C For each dimension, we need to calculate powers of an
C appropriate irreducible polynomial :  see Niederreiter
C page 65, just below equation (19).
C   Read the appropriate irreducible polynomial into PX,
C and its degree into E.  Set polynomial B = PX ** 0 = 1.
C M is the degree of B.  Subsequently B will hold higher
C powers of PX.
C   The polynomial PX is stored in file 'irrtabs.txt' in the
C format
C   n  a0  a1  a2  ... an
C where n is the degree of the polynomial and the ai are
C its coefficients.
C
        READ (GFUNIT, 900) E, (PX(J), J = 0,E)
        PX(DEG) = E
        B(DEG) = 0
        B(0) = 1
C
C Niederreiter (page 56, after equation (7), defines two
C variables Q and U.  We do not need Q explicitly, but we
C do need U.
C
        U = 0
C
        DO 90 J = 1, NFIGS
C
C  If U = 0, we need to set B to the next power of PX
C  and recalculate V.  This is done by subroutine CALCV.
C
          IF (U .EQ. 0) CALL CALCV (PX, B, V, MAXV)
C
C Now C is obtained from V.  Neiderreiter
C obtains A from V (page 65, near the bottom), and
C then gets C from A (page 56, equation (7)).
C However this can be done in one step.
C
          DO 50 R = 0, NFIGS-1
            C(I,J,R) = V(R+U)
   50     CONTINUE
C
C Increment U.  If U = E, then U = 0 and in Niederreiter's
C paper Q = Q + 1.  Here, however, Q is not used explicitly.
C
          U = U + 1
          IF (U .EQ. E) U = 0
  90    CONTINUE
C
 1000 CONTINUE
C
      CLOSE (GFUNIT)
      RETURN
C
  500 WRITE (*,*) ' CALCC :  Tables for q =', Q, ' not found'
      STOP
      END
C
C     ***** end of SUBROUTINE CALCC
      SUBROUTINE CALCV (PX, B, V, MAXV)
C
C   This version :  12 February 1991
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This program calculates the values of the constants V(J,R) as
C   described in BFN section 3.3.  It is called from either CALCC or
C   CALCC2.  The values transmitted through common /FIELD/ determine
C   which field we are working in.
C
C INPUT :
C   PX is the appropriate irreducible polynomial for the dimension
C     currently being considered.  The degree of PX will be called E.
C   B is the polynomial defined in section 2.3 of BFN.  On entry to
C     the subroutine, the degree of B implicitly defines the parameter
C     J of section 3.3, by degree(B) = E*(J-1).
C   MAXV gives the dimension of the array V.
C   On entry, we suppose that the common block /FIELD/ has been set
C     up correctly (using SETFLD).
C
C OUTPUT :
C   On output B has been multiplied by PX, so its degree is now E*J.
C   V contains the values required.
C
C USES :
C   The subroutine PLYMUL is used to multiply polynomials.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication, and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER MAXV, E, I, J, KJ, M, BIGM, R, TERM
      INTEGER PX(-1:MAXDEG), B(-1:MAXDEG), V(0:MAXV)
      INTEGER H(-1:MAXDEG)
C
      INTEGER ARBIT, NONZER
      ARBIT() = 1
C
C   We use ARBIT() to indicate where the user can place
C   an arbitrary element of the field of order Q, while NONZER
C   shows where he must put an arbitrary non-zero element
C   of the same field.  For the code,
C   this means 0 <= ARBIT < Q and 0 < NONZER < Q.  Within these
C   limits, the user can do what he likes.  ARBIT is declared as
C   a function as a reminder that a different arbitrary value may
C   be returned each time ARBIT is referenced.
C
C    BIGM is the M used in section 3.3.
C    It differs from the [little] m used in section 2.3,
C    denoted here by M.
C
      NONZER = 1
C
      E = PX(DEG)
C
C   The poly H is PX**(J-1), which is the value of B on arrival.
C   In section 3.3, the values of Hi are defined with a minus sign :
C   don't forget this if you use them later !
C
      DO 10 I = -1, B(DEG)
   10   H(I) = B(I)
      BIGM = H(DEG)
C
C   Now multiply B by PX so B becomes PX**J.
C   In section 2.3, the values of Bi are defined with a minus sign :
C   don't forget this if you use them later !
C
      CALL PLYMUL (PX, B, B)
      M = B(DEG)
C
C   We don't use J explicitly anywhere, but here it is just in case.
C
      J = M / E
C
C   Now choose a value of Kj as defined in section 3.3.
C   We must have 0 <= Kj < E*J = M.
C   The limit condition on Kj does not seem very relevant
C   in this program.
C
      KJ = BIGM
C
C   Now choose values of V in accordance with the conditions in
C   section 3.3
C
      DO 20 R = 0, KJ-1
   20   V(R) = 0
      V(KJ) = 1
C
      IF (KJ .LT. BIGM) THEN
C
        TERM = SUB (0, H(KJ))
C
        DO 30 R = KJ+1, BIGM-1
          V(R) = ARBIT()
C
C         Check the condition of section 3.3,
C         remembering that the H's have the opposite sign.
C
          TERM = SUB (TERM, MUL (H(R), V(R)))
   30   CONTINUE
C
C         Now V(BIGM) is anything but TERM
C
          V(BIGM) = ADD (NONZER, TERM)
C
        DO 40 R = BIGM+1, M-1
   40     V(R) = ARBIT()
C
      ELSE
C       This is the case KJ .GE. BIGM
C
        DO 50 R = KJ+1, M-1
   50     V(R) = ARBIT()
C
      ENDIF
C
C   Calculate the remaining V's using the recursion of section 2.3,
C   remembering that the B's have the opposite sign.
C
      DO 70 R = 0, MAXV-M
        TERM = 0
        DO 60 I = 0, M-1
          TERM = SUB (TERM, MUL (B(I), V(R+I)))
   60   CONTINUE
        V(R+M) = TERM
   70 CONTINUE
C
      RETURN
      END
C
C     ***** end of SUBROUTINE CALCV
      INTEGER FUNCTION CHARAC (QIN)
C
C   This version :  12 December 1991
C
C   This function gives the characteristic for a field of
C   order QIN.  If no such field exists, or if QIN is out of
C   the range we can handle, returns 0.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER QIN, CH(MAXQ)
      SAVE CH
C
      DATA CH  / 0,  2,  3,  2,  5,  0,  7,  2,  3,  0,
     1          11,  0, 13,  0,  0,  2, 17,  0, 19,  0,
     2           0,  0, 23,  0,  5,  0,  3,  0, 29,  0,
     3          31,  2,  0,  0,  0,  0, 37,  0,  0,  0,
     4          41,  0, 43,  0,  0,  0, 47,  0,  7,  0/
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        CHARAC = 0
      ELSE
        CHARAC = CH(QIN)
      ENDIF
C
      END
C
C     ***** end of INTEGER FUNCTION CHARAC
      SUBROUTINE SETFLD (QIN)
      INTEGER QIN
C
C   This version : 12 December 1991
C
C   This subroutine sets up addition, multiplication, and
C   subtraction tables for the finite field of order QIN.
C   If necessary, it reads precalculated tables from the file
C   'gftabs.txt' using unit 1.  These precalculated tables
C   are supposed to have been put there by GFARIT.
C
C      *****  For the base-2 programs, these precalculated
C      *****  tables are not needed and, therefore, neither
C      *****  is GFARIT.
C
C
C   Unit 1 is closed both before and after the call of SETFLD.
C
C USES
C   Integer function CHARAC gets the characteristic of a field.
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, N, CHARAC
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' SETFLD :  Bad value of Q'
        STOP
      ENDIF
C
      Q = QIN
      P = CHARAC(Q)
C
      IF (P .EQ. 0) THEN
        WRITE (*,*) ' SETFLD :  There is no field of order', Q
        STOP
      ENDIF
C
C Set up to handle a field of prime order :  calculate ADD and MUL.
C
      IF (P .EQ. Q) THEN
        DO 10 I = 0, Q-1
          DO 10 J = 0, Q-1
            ADD(I,J) = MOD(I+J, P)
            MUL(I,J) = MOD(I*J, P)
   10   CONTINUE
C
C Set up to handle a field of prime-power order :  tables for
C ADD and MUL are in the file 'gftabs.txt'.
C
      ELSE
        OPEN (UNIT=1, FILE='gftabs.txt', STATUS='old')
C
C    *****  OPEN statements are system dependent.
C
   20   READ (1, 900, END=500) N
  900   FORMAT (20I3)
        DO 30 I = 0, N-1
          READ (1, 900) (ADD(I,J), J = 0, N-1)
   30   CONTINUE
        DO 40 I = 0, N-1
          READ (1, 900) (MUL(I,J), J = 0, N-1)
   40   CONTINUE
        IF (N .NE. Q) GOTO 20
        CLOSE (1)
      ENDIF
C
C Now use the addition table to set the subtraction table.
C
      DO 60 I = 0, Q-1
        DO 50 J = 0, Q-1
          SUB(ADD(I,J), I) = J
   50   CONTINUE
   60 CONTINUE
      RETURN
C
  500 WRITE (*,*) ' SETFLD :  Tables for q =', Q, ' not found'
      STOP
C
      END
C
C     ***** end of SUBROUTINE SETFLD
      SUBROUTINE PLYMUL (PA, PB, PC)
C
C   This version :  12 December 1991
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, DEGA, DEGB, DEGC, TERM
      INTEGER PA(-1:MAXDEG), PB(-1:MAXDEG), PC(-1:MAXDEG)
      INTEGER PT(-1:MAXDEG)
C
C   Multiplies polynomial PA by PB putting the result in PC.
C   Coefficients are elements of the field of order Q.
C
      DEGA = PA(DEG)
      DEGB = PB(DEG)
      IF (DEGA .EQ. -1 .OR. DEGB .EQ. -1) THEN
        DEGC = -1
      ELSE
        DEGC = DEGA + DEGB
      ENDIF
      IF (DEGC .GT. MAXDEG) THEN
        WRITE (*,*) ' PLYMUL :  Degree of product exceeds MAXDEG'
        STOP
      ENDIF
C
      DO 20 I = 0, DEGC
        TERM = 0
        DO 10 J = MAX(0, I-DEGA), MIN(DEGB, I)
   10     TERM = ADD(TERM, MUL(PA(I-J), PB(J)))
   20   PT(I) = TERM
C
      PC(DEG) = DEGC
      DO 30 I = 0, DEGC
   30   PC(I) = PT(I)
      DO 40 I = DEGC+1, MAXDEG
   40   PC(I) = 0
      RETURN
      END
C
C     *****   end of SUBROUTINE PLYMUL
      SUBROUTINE GOLO (QUASI)
      REAL QUASI(*)
C
C This version : 21 February 1992
C
C See the general comments on implementing Niederreiter's
C low-discrepancy sequences.
C
C Call subroutine INLO before calling GOLO.  Thereafter
C GOLO generates a new quasi-random vector on each call,
C
C INPUT
C   From INLO, labelled common /COMM/ and labelled common
C     /FIELD/, properly initialized.  (INLO calls SETFLD
C     to initialize /FIELD/).
C
C OUTPUT
C   To the user's program, the next vector in the sequence in
C     array QUASI.
C
C   -------------------------------------------------------------
C
C
C   This segment defines the common block /COMM/ and some
C   associated parameters.  These are for use in the general base
C   version of the generator.
C
      INTEGER MAXDIM, MAXFIG, NBITS
      PARAMETER (MAXDIM=12, MAXFIG=20, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   MAXFIG is the maximum number of base-Q digits we can handle.
C   MAXINT is the largest fixed point integer we can represent.
C   NBITS is the number of bits in a fixed-point integer, not
C     counting the sign.
C     ***** NBITS is machine dependent *****
C
      INTEGER C(MAXDIM, MAXFIG, 0:MAXFIG-1)
      INTEGER COUNT(0:MAXFIG-1), D(MAXDIM, MAXFIG)
      INTEGER NEXTQ(MAXDIM), QPOW(MAXFIG)
      INTEGER DIMEN, NFIGS
      REAL    RECIP
      COMMON  /COMM/ C, COUNT, D, NEXTQ, QPOW, DIMEN, NFIGS, RECIP
      SAVE    /COMM/
C
C   The common block /COMM/ :
C     C     - Contains the values of Niederreiter's C(I,J,R)
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of base-Q digits.  COUNT(R)
C             is the same as Niederreiter's AR(N) (page 54)
C             except that N is implicit.
C     D     - The values of D(I,J).
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP.
C     QPOW  - To speed things up a bit. QPOW(I) = Q ** (NFIGS-I).
C     DIMEN   - The dimension of the sequence to be generated
C     NFIGS - The number of base-Q digits we are using.
C     RECIP - 1.0 / (Q ** NFIGS)
C
C   Array C of the common block is set up by subroutine CALCC.
C   The other items in the common block are set up by INLO.
C
C   ------------------------------------------------------------
C
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, R, OLDCNT, DIFF, NQ
C
C Multiply the numerators in NEXTQ by RECIP to get the next
C quasi-random vector
C
      DO 5 I = 1, DIMEN
        QUASI(I) = NEXTQ(I) * RECIP
    5 CONTINUE
C
C Update COUNT, treated as a base-Q integer.  Instead of
C recalculating the values of D from scratch, we update
C them for each digit of COUNT which changes.  In terms of
C Niederreiter page 54, NEXTQ(I) corresponds to XI(N), with
C N being implicit, and D(I,J) corresponds to XI(N,J), again
C with N implicit.  Finally COUNT(R) corresponds to AR(N).
C
      R = 0
   10 IF (R .EQ. NFIGS) THEN
        WRITE (*,*) ' Too many calls on subroutine GOLO'
        STOP
      ENDIF
      OLDCNT = COUNT(R)
      IF (COUNT(R) .LT. Q-1) THEN
        COUNT(R) = COUNT(R) + 1
      ELSE
        COUNT(R) = 0
      ENDIF
      DIFF = SUB(COUNT(R), OLDCNT)
C
C Digit R has just changed.  DIFF says how much it changed
C by.  We use this to update the values of array D.
C
      DO 40 I = 1, DIMEN
        DO 30 J = 1, NFIGS
          D(I,J) = ADD(D(I,J), MUL(C(I,J,R), DIFF))
   30   CONTINUE
   40 CONTINUE
C
C If COUNT(R) is now zero, we must propagate the carry
C
      IF (COUNT(R).EQ.0) THEN
        R = R + 1
        GOTO 10
      ENDIF
C
C Now use the updated values of D to calculate NEXTQ.
C Array QPOW helps to speed things up a little :
C   QPOW(J) is Q ** (NFIGS-J).
C
      DO 60 I = 1, DIMEN
        NQ = 0
        DO 50 J = 1, NFIGS
          NQ = NQ + D(I,J) * QPOW(J)
   50   CONTINUE
        NEXTQ(I) = NQ
   60 CONTINUE
C
      RETURN
      END
C
C     ***** end of SUBROUTINE GOLO
      REAL FUNCTION TESTF (N, DIMEN, QUASI)
      INTEGER I, N, DIMEN
      REAL X, EXACTF, QUASI(*)
C
C This version :  4 Mar 1992
C
C Provides a variety of test integrals for quasi-random
C sequences.  A call on TESTF computes an estimate of the
C integral ;  a call on EXACTF computes the exact value.
C
      GOTO (100, 200, 300, 400) N
C
      ENTRY EXACTF (N, DIMEN)
      GOTO (1100, 1200, 1300, 1400) N
C
C Test integral 1
C
  100 TESTF = 1.0
      DO 110 I = 1, DIMEN
        TESTF = TESTF * ABS(4 * QUASI(I) - 2)
  110 CONTINUE
      RETURN
C
 1100 EXACTF = 1.0
      RETURN
C
C Test integral 2
C
  200 TESTF = 1.0
      DO 210 I = 1, DIMEN
        TESTF = TESTF * I * COS(I * QUASI(I))
  210 CONTINUE
      RETURN
C
 1200 EXACTF = 1.0
      DO 1210 I = 1, DIMEN
        EXACTF = EXACTF * SIN(FLOAT(I))
 1210 CONTINUE
      RETURN
C
C Test integral 3
C
  300 TESTF = 1.0
      DO 350 I = 1, DIMEN
        X = 2 * QUASI(I) - 1
        GOTO (310, 320, 330, 340) MOD(I, 4)
  310   TESTF = TESTF * X
        GOTO 350
  320   TESTF = TESTF * (2*X*X - 1)
        GOTO 350
  330   TESTF = TESTF * (4*X*X - 3) * X
        GOTO 350
  340   X = X * X
        TESTF = TESTF * (8*X*X - 8*X + 1)
  350 CONTINUE
      RETURN
C
 1300 EXACTF = 0.0
      RETURN
C
C Test integral 4
C
  400        TESTF = 0
             X = 1
             DO 410   I = 1, DIMEN
                  X = - X * QUASI(I)
                  TESTF = TESTF + X
  410        CONTINUE
             RETURN
C
C
 1400 X = 1.0 / (2 ** (DIMEN ))
      IF (MOD(DIMEN, 2) .EQ. 0) THEN
        EXACTF = (X - 1) / 3
      ELSE
        EXACTF = (X + 1) / 3
      ENDIF
      RETURN
C
      END
