c*** gfplys.f
      PROGRAM GFPLYS
C
C This version :  12 December 1991
C
C   The program calculates irreducible polynomials for various
C   finite fields, and writes them out to the file 'irrtabs.txt'.
C   Finite field arithmetic is carried out with the help of
C   precalculated addition and multiplication tables found on
C   the file 'gftabs.txt'.  The format of the irreducible polys on
C   the output file is
C       Q
C       d1   a1  a2 ... a(d1)
C       d2   b1  b2 ... b(d2)
C           etc
C   where Q is the order of the field, d1 is the degree of the
C   first irreducible poly, a1, a2, ..., a(d1) are its
C   coefficients, and so on.
C
C   The file 'gftabs.txt' is read on unit 1.  GFPLYS and the
C   associated subroutines assume that GFARIT has been run previously
C   to put the required data in this file.  GFPLYS writes its
C   output on file 'irrtabs.txt' using unit 2.
C
C   GFPLYS and its associated subroutines are listed below.
C   An asterisk indicates a subroutine also used by GFARIT
C   and by GENIN.  An ampersand indicates a subroutine also
C   used by GFARIT but not by GENIN.
C
C      GFPLYS
C         IRRED       (reads unit 1, writes unit 2)
C         CHARAC *
C         SETFLD *
C         ITOP   &
C         PTOI   &
C         PLYMUL *
C         FIND   
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
      INTEGER I
C
      OPEN (UNIT=1, FILE='gftabs.txt', STATUS = 'old')
      OPEN (UNIT=2, FILE='irrtabs.txt', STATUS = 'unknown')
C
C      *****  OPEN statements are system dependent
C
      WRITE (*,*) ' This is GFPLYS'
      DO 10 I = 2, MAXQ
        CALL IRRED(I)
        WRITE (*,*) ' GFPLYS :  Case ', I, ' done'
   10 CONTINUE
C
      END
C
C     *****  end of PROGRAM GFPLYS
      SUBROUTINE IRRED (QIN)
C
C   This version :  12 December 1991
C
C
C   This routine reads the file gftabs.txt from unit 1
C   and writes the file irrtabs.txt onto unit 2.
C   GFPLYS opens units 1 and 2.
C   SETFLD opens unit 1, reads gftabs.txt, and then closes unit 1
C   on each call.
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
      INTEGER NPOLS, MAXS
      PARAMETER (NPOLS=25, MAXS=400)
      LOGICAL SEIVE(MAXS)
      INTEGER MONPOL(MAXS)
C
C   The parameter NPOLS says how many irreducible polys are to
C   be calculated for each field.
C   We find the irreducible polys using a sieve.  Parameter
C   MAXS gives the size of this seive.  Array MONPOL holds monic
C   polys, array SEIVE says whether the poly is still OK.
C
      INTEGER QIN, I, J, K, N, PTOI, FIND, CHARAC
      INTEGER PI(-1:MAXDEG), PJ(-1:MAXDEG), PK(-1:MAXDEG)
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' IRRED :  Bad value of Q'
        RETURN
      ENDIF
C
      P = CHARAC(QIN)
C
C   If no field of order QIN exists, there is nothing to do.
C
      IF (P .EQ. 0) RETURN
C
C   Otherwise, set up the field arithmetic tables.
C
      CALL SETFLD (QIN)
C
C   Set up seive containing only monic polys
C
      I = 0
      J = 1
      K = Q
      DO 50 N = 1, MAXS
        I = I + 1
        IF (I .EQ. J) THEN
          I = K
          J = 2 * K
          K = Q * K
        ENDIF
        MONPOL(N) = I
        SEIVE(N) = .TRUE.
   50 CONTINUE
C
C   Write out irreducible polys as they are found
C
      N = 0
      WRITE (2, 900) QIN
  900 FORMAT (20I3)
      DO 200 I = 1, MAXS
        IF (SEIVE(I)) THEN
          CALL ITOP (MONPOL(I), PI)
          K = PI(DEG)
          WRITE (2, 900) K, (PI(J), J = 0, K)
          N = N + 1
          IF (N .EQ. NPOLS) RETURN
C
          DO 100 J = I, MAXS
            CALL ITOP (MONPOL(J), PJ)
            CALL PLYMUL (PI, PJ, PK)
            K = FIND (PTOI (PK), MONPOL, J, MAXS)
            IF (K .NE. 0) SEIVE(K) = .FALSE.
  100     CONTINUE
        ENDIF
  200 CONTINUE
C
      WRITE (*,*) ' IRRED :  Seive too small'
      WRITE (*,*) ' Only', N, ' irreducible polys were found'
      RETURN
C
      END
C
C     *****   end of SUBROUTINE IRRED
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
      SUBROUTINE ITOP (IN, POLY)
C
C   This version :  12 December 1991
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
      INTEGER IN, I, J, POLY(-1:MAXDEG)
C
C   Converts an integer to a polynomial with coefficients in the
C   field of order Q.
C
      DO 10 J = -1, MAXDEG
   10   POLY(J) = 0
C
      I = IN
      J = -1
   20 IF (I .GT. 0) THEN
        J = J + 1
        IF (J .GT. MAXDEG) THEN
          WRITE (*,*) ' ITOP :  Polynomial exceeds MAXDEG'
          STOP
        ENDIF
        POLY(J) = MOD (I, Q)
        I = I / Q
        GOTO 20
      ENDIF
      POLY(DEG) = J
      RETURN
      END
C
C     *****  end of SUBROUTINE ITOP
      INTEGER FUNCTION PTOI (POLY)
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
      INTEGER I, J, POLY(-1:MAXDEG)
C
C   Converts a polynomial with coefficients in the field of
C   order Q to an integer.
C
      I = 0
      DO 10 J = POLY(DEG), 0, -1
   10   I = Q * I + POLY(J)
      PTOI = I
      RETURN
      END
C
C     *****  end of INTEGER FUNCTION PTOI
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
      INTEGER FUNCTION FIND (N, TAB, I, MAXTAB)
C
C   This version :  12 December 1991
C
      INTEGER N, I, MAXTAB, TAB(MAXTAB), J
C
C   Look up N in ordered TAB(I) to TAB(MAXTAB)
C
      FIND = 0
      IF (N .GT. TAB(MAXTAB)) RETURN
      DO 10 J = I, MAXTAB
        IF (TAB(J) .EQ. N) THEN
          FIND = J
          RETURN
        ENDIF
   10 CONTINUE
      RETURN
      END
