c*** gfarit.f
      PROGRAM GFARIT
C
C   This version :  12 December 1991
C
C   The program calculates addition and multiplication tables
C   for arithmetic in finite fields, and writes them out to
C   the file 'gftabs.txt'.  Tables are only calculated for fields
C   of prime-power order Q, the other cases being trivial.
C   For each value of Q, the file contains first Q, then the
C   addition table, and lastly the multiplication table.
C
C
C   GFARIT and its associated subroutines below must be
C   run to set up the file gftabs.txt [on unit 1].  
C   After gftabs.txt has been set up, run GFPLYS with its
C   associated subroutine to set up the file irrtabs.txt
C   [on unit 2]; this requires
C   reading gftabs.txt from unit 1.  The files gftabs.txt
C   and irrtabs.txt can be saved for future use.  Thus, each
C   user [or set of users with access to these files] needs to
C   run the respective sets of programs associated with GFARIT
C   and GFPLYS just once.  This must be done before running the
C   set of programs associated with GENIN.  The set of programs
C   tailored for base 2, using the driver GENIN2, requires
C   neither gftabs.txt nor irrtabs.txt, hence neither GFARIT 
C   nor GFPLYS.
C
C   Below we list [the main program] GFARIT and its associated
C   subroutines.  An asterisk indicates a subroutine also used
C   by GFPLYS and by GENIN.  An ampersand denotes a subroutine
C   also used by GFPLYS but not by GENIN.
C
C       GFARIT
C          GFTAB           (writes unit 1)
C          CHARAC *
C          SETFLD *
C          ITOP   &
C          PTOI   &
C          PLYADD 
C          PLYMUL *
C          PLYDIV
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
      INTEGER I
C
      OPEN (UNIT=1, FILE='gftabs.txt', STATUS = 'unknown')
C
C       ******  OPEN statements are system dependent
C
      WRITE (*,*) ' This is GFARIT'
      DO 10 I = 2, MAXQ
        CALL GFTAB(I)
        WRITE (*,*) ' GFARIT :  Case ', I, ' done'
   10 CONTINUE
C
      WRITE (*,*) ' GFARIT ends'
      STOP
      END
C
C     *****  end of PROGRAM GFARIT
      SUBROUTINE GFTAB (QIN)
C
C   This version :  12 December 1991
C
C
C    This routine writes gftabs.txt onto unit 1.
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
C   The common /FIELD/ will be set up to work modulo P, the
C   characteristic of field QIN.  We construct the tables for
C   the field of order QIN in GFADD and GFMUL.
C
      INTEGER QIN, I, J, PTOI, CHARAC
      INTEGER PI(-1:MAXDEG), PJ(-1:MAXDEG), PK(-1:MAXDEG)
      INTEGER GFADD(0:MAXQ-1, 0:MAXQ-1), GFMUL(0:MAXQ-1, 0:MAXQ-1)
C
C   IRRPLY holds irreducible polynomials for constructing
C   prime-power fields.  IRRPLY(-2,I) says which field this
C   row is used for, and then the rest of the row is a
C   polynomial (with the degree in IRRPLY(-1,I) as usual).
C   The chosen irreducible poly is copied into MODPLY for use.
C
      INTEGER IRRPLY(8, -2:MAXDEG), MODPLY(-1:MAXDEG)
      SAVE IRRPLY
C
      DATA (IRRPLY(1,J), J=-2,2) /4, 2, 1, 1, 1/
      DATA (IRRPLY(2,J), J=-2,3) /8, 3, 1, 1, 0, 1/
      DATA (IRRPLY(3,J), J=-2,2) /9, 2, 1, 0, 1/
      DATA (IRRPLY(4,J), J=-2,4) /16, 4, 1, 1, 0, 0, 1/
      DATA (IRRPLY(5,J), J=-2,2) /25, 2, 2, 0, 1/
      DATA (IRRPLY(6,J), J=-2,3) /27, 3, 1, 2, 0, 1/
      DATA (IRRPLY(7,J), J=-2,5) /32, 5, 1, 0, 1, 0, 0, 1/
      DATA (IRRPLY(8,J), J=-2,2) /49, 2, 1, 0, 1/
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' ARITH :  Bad value of Q'
        RETURN
      ENDIF
C
      P = CHARAC(QIN)
C
C   If QIN is not a prime-power, we are not interested.
C
      IF (P .EQ. 0 .OR. P .EQ. QIN) RETURN
C
C   Otherwise, we set up the elements of the common /FIELD/
C   ready to do arithmetic mod P, the characteristic of QIN.
C
      CALL SETFLD ( QIN )
C
C   Next find a suitable irreducible polynomial and copy it
C   to array MODPLY.
C
      I = 1
   20 IF (IRRPLY(I,-2) .NE. QIN) THEN
        I = I + 1
        GOTO 20
      ENDIF
      DO 30 J = -1, IRRPLY(I, DEG)
   30   MODPLY(J) = IRRPLY(I, J)
      DO 40 J = IRRPLY(I,DEG)+1, MAXDEG
   40   MODPLY(J) = 0
C
C   Deal with the trivial cases ...
C
      DO 60 I = 0, QIN-1
        GFADD(I,0) = I
        GFADD(0,I) = I
        GFMUL(I,0) = 0
   60   GFMUL(0,I) = 0
      DO 70 I = 1, QIN-1
        GFMUL(I,1) = I
   70   GFMUL(1,I) = I
C
C   ... and now deal with the rest.  Each integer from 1 to QIN-1
C   is treated as a polynomial with coefficients handled mod P.
C   Multiplication of polynomials is mod MODPLY.
C
      DO 80 I = 1, QIN-1
        CALL ITOP (I, PI)
        DO 80 J = 1, I
          CALL ITOP (J, PJ)
          CALL PLYADD (PI, PJ, PK)
          GFADD(I,J) = PTOI (PK)
          GFADD(J,I) = GFADD(I,J)
          IF (I .GT. 1 .AND. J .GT. 1) THEN
            CALL PLYMUL (PI, PJ, PK)
            CALL PLYDIV (PK, MODPLY, PJ, PK)
            GFMUL(I,J) = PTOI (PK)
            GFMUL(J,I) = GFMUL(I,J)
          ENDIF
   80 CONTINUE
C
C   Write the newly-calculated tables out to unit 1
C
      WRITE (1, 900) QIN
C
C   This is the table gftabs.txt.
C   SETFLD opens unit 1, reads gftabs.txt, and then closes unit 1.
C
      DO 100 I = 0, QIN-1
  100   WRITE (1, 900) (GFADD(I,J), J = 0, QIN-1)
      DO 110 I = 0, QIN-1
  110   WRITE (1, 900) (GFMUL(I,J), J = 0, QIN-1)
  900 FORMAT (20I3)
C
      RETURN
      END
C
C     *****  end of SUBROUTINE GFTAB
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
      DO 10 I = 0, Q-1
        DO 10 J = 0, Q-1
          ADD(I,J) = MOD(I+J, P)
          MUL(I,J) = MOD(I*J, P)
   10 CONTINUE
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
      SUBROUTINE PLYADD (PA, PB, PC)
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
      INTEGER I, MAXAB, DEGC
      INTEGER PA(-1:MAXDEG), PB(-1:MAXDEG), PC(-1:MAXDEG)
C
C   Adds polynomials PA and PB putting result in PC.
C   Coefficients are elements of the field of order Q.
C
      MAXAB = MAX (PA(DEG), PB(DEG))
      DEGC = -1
      DO 10 I = 0, MAXAB
        PC(I) = ADD (PA(I), PB(I))
        IF (PC(I) .NE. 0) DEGC = I
   10 CONTINUE
      PC(DEG) = DEGC
      DO 20 I = MAXAB+1, MAXDEG
   20   PC(I) = 0
      RETURN
      END
C
C     *****   end of SUBROUTINE PLYADD
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
      SUBROUTINE PLYDIV (PA, PB, PQ, PR)
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
      INTEGER I, J, D, M, BINV, DEGB, DEGR, DEGQ
      INTEGER PA(-1:MAXDEG), PB(-1:MAXDEG), PQ(-1:MAXDEG)
      INTEGER PR(-1:MAXDEG)
      INTEGER PTQ(-1:MAXDEG), PTR(-1:MAXDEG)
C
C   Divides polynomial PA by PB, putting the quotient in PQ
C   and the remainder in PR.
C   Coefficients are elements of the field of order Q.
C
      IF (PB(DEG) .EQ. -1) THEN
        WRITE (*,*) ' PLYDIV :  Divide by zero'
        STOP
      ENDIF
C
      DO 10 I = -1, MAXDEG
        PTQ(I) = 0
   10   PTR(I) = PA(I)
      DEGR = PA(DEG)
      DEGB = PB(DEG)
      DEGQ = DEGR - DEGB
      IF (DEGQ .LT. 0) DEGQ = -1
C
C   Find the inverse of the leading coefficient of PB.
C
      J = PB(DEGB)
      DO 15 I = 1, P-1
        IF (MUL(I,J) .EQ. 1) BINV = I
   15 CONTINUE
C
      DO 30 D = DEGQ, 0, -1
        M = MUL (PTR(DEGR), BINV)
        DO 20 I = DEGB, 0, -1
   20     PTR(DEGR+I-DEGB) = SUB (PTR(DEGR+I-DEGB), MUL (M, PB(I)))
        DEGR = DEGR - 1
   30   PTQ(D) = M
C
      DO 40 I = 0, MAXDEG
        PQ(I) = PTQ(I)
   40   PR(I) = PTR(I)
      PQ(DEG) = DEGQ
   50 IF (PR(DEGR) .EQ. 0 .AND. DEGR .GE. 0) THEN
        DEGR = DEGR - 1
        GOTO 50
      ENDIF
      PR(DEG) = DEGR
C
      RETURN
      END
