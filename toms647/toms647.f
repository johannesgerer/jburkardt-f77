      block data bdhalt

c*********************************************************************72
c
cc BDHALT initializes labeled common for INHALT and GOHALT.
c
      integer s
      double precision e
      double precision prime(40)

      common /halton/ e, prime, s

      save /halton/

      data (prime(i),i=1,40) / 
     &    2.0D+00,   3.0D+00,   5.0D+00,   7.0D+00,  11.0D+00,
     &   13.0D+00,  17.0D+00,  19.0D+00,  23.0D+00,  29.0D+00,
     &   31.0D+00,  37.0D+00,  41.0D+00,  43.0D+00,  47.0D+00,
     &   53.0D+00,  59.0D+00,  61.0D+00,  67.0D+00,  71.0D+00,
     &   73.0D+00,  79.0D+00,  83.0D+00,  89.0D+00,  97.0D+00,
     &  101.0D+00, 103.0D+00, 107.0D+00, 109.0D+00, 113.0D+00,
     &  127.0D+00, 131.0D+00, 137.0D+00, 139.0D+00, 149.0D+00,
     &  151.0D+00, 157.0D+00, 163.0D+00, 167.0D+00, 173.0D+00 /

      end
      block data bdsobl

c*********************************************************************72
C
cc BDSOBL initializes labeled common for INSOBL.
C
C     THE ARRAY POLY GIVES SUCCESSIVE PRIMITIVE
C     POLYNOMIALS CODED IN BINARY, E.G.
C          45 = 100101
C     HAS BITS 5, 2, AND 0 SET (COUNTING FROM THE
C     RIGHT) AND THEREFORE REPRESENTS
C          X**5 + X**2 + X**0
C
C     THESE  POLYNOMIALS ARE IN THE ORDER USED BY
C     SOBOL IN USSR COMPUT. MATHS. MATH. PHYS. 16 (1977),
C     236-242. A MORE COMPLETE TABLE IS GIVEN IN SOBOL AND
C     LEVITAN, THE PRODUCTION OF POINTS UNIFORMLY
C     DISTRIBUTED IN A MULTIDIMENSIONAL CUBE (IN RUSSIAN),
C     PREPRINT IPM AKAD. NAUK SSSR, NO. 40, MOSCOW 1976.
C
C         THE INITIALIZATION OF THE ARRAY V IS FROM THE
C     LATTER PAPER. FOR A POLYNOMIAL OF DEGREE M, M INITIAL
C     VALUES ARE NEEDED :  THESE ARE THE VALUES GIVEN HERE.
C     SUBSEQUENT VALUES OF V ARE CALCULATED IN "INSOBL".
C
      INTEGER V(40,30),S,MAXCOL,COUNT
      integer lastq(40)
      integer POLY(40),I
      REAL RECIPD

      COMMON /SOBOL/ V, S, MAXCOL, COUNT, POLY, LASTQ, RECIPD

      SAVE /SOBOL/

      DATA (POLY(I),I=1,40)/1,3,7,11,13,19,25,37,59,47,61,55,41,67,97,
     .     91,109,103,115,131,193,137,145,143,241,157,185,167,229,171,
     .     213,191,253,203,211,239,247,285,369,299/

      DATA (V(I,1),I=1,40)/40*1/
      DATA (V(I,2),I=3,40)/1,3,1,3,1,3,3,1,3,1,3,1,3,1,1,3,1,3,1,3,1,3,
     .     3,1,3,1,3,1,3,1,1,3,1,3,1,3,1,3/
      DATA (V(I,3),I=4,40)/7,5,1,3,3,7,5,5,7,7,1,3,3,7,5,1,1,5,3,3,1,7,
     .     5,1,3,3,7,5,1,1,5,7,7,5,1,3,3/
      DATA (V(I,4),I=6,40)/1,7,9,13,11,1,3,7,9,5,13,13,11,3,15,5,3,15,7,
     .     9,13,9,1,11,7,5,15,1,15,11,5,3,1,7,9/
      DATA (V(I,5),I=8,40)/9,3,27,15,29,21,23,19,11,25,7,13,17,1,25,29,
     .     3,31,11,5,23,27,19,21,5,1,17,13,7,15,9,31,9/
      DATA (V(I,6),I=14,40)/37,33,7,5,11,39,63,27,17,15,23,29,3,21,13,
     .     31,25,9,49,33,19,29,11,19,27,15,25/
      DATA (V(I,7),I=20,40)/13,33,115,41,79,17,29,119,75,73,105,7,59,65,
     .     21,3,113,61,89,45,107/
      DATA (V(I,8),I=38,40)/7,23,39/

      end
      function exor ( iin, jin )

c*********************************************************************72
C
cc EXOR calculates the exclusive OR of its input arguments.
c
C     THIS FUNCTION CALCULATES THE EXCLUSIVE-OR OF ITS
C     TWO INPUT PARAMETERS
C
      integer exor
      INTEGER I,J,K,L,IIN,JIN,I2,J2

      I = IIN
      J = JIN
      K = 0
      L = 1

   10 continue

      IF (I.EQ.0 .AND. J.EQ.0) THEN
          EXOR = K
          RETURN
      END IF
C
C     CHECK THE CURRENT RIGHT-HAND BITS OF I AND J.
C     IF THEY DIFFER, SET THE APPROPRIATE BIT OF K.
C
      I2 = I/2
      J2 = J/2
      IF ((I.EQ.2*I2) .NEQV. (J.EQ.2*J2)) K = K + L
      I = I2
      J = J2
      L = 2*L
      GO TO 10

      end
      subroutine gofaur ( quasi )

c*********************************************************************72
c
cc GOFAUR generates a new Faure quasirandom vector with each call.
c
C       IT IMPLEMENTS A METHOD OF H.FAURE,
C       "ACTA ARITHMETICA XLI(1982),337-351".
C       (SEE ESPECIALLY PAGE 342).
C
C       THE USER MUST CALL "INFAUR" BEFORE
C       CALLING "GOFAUR".
C       AFTER CALLING "INFAUR", TEST FLAG(1)
C       AND FLAG(2); IF EITHER IS FALSE, DO
C       NOT CALL GOFAUR. READ THE COMMENTS AT
C       THE BEGINNING OF INFAUR AND THEN
C       THOSE BELOW.
C
C       ALL INPUTS COME FROM "INFAUR" VIA
C       LABELLED COMMON "FAURE"; FOR THEIR
C       DEFINITIONS, SEE "INFAUR".
C
C       INPUTS:
C         S,QS,COEF,NEXTN,TESTN,HISUM,RQS
C
C       OUTPUTS:
C         TO USER'S CALLING PROGRAM:
C         QUASI - A NEW QUASIRANDOM VECTOR
C
      INTEGER S,QS,COEF(0:19,0:19),NEXTN,TESTN,HISUM,I,J,K,YTEMP(0:19),
     .        ZTEMP,KTEMP,LTEMP,MTEMP
C
      REAL QUASI(40),RQS,R
C
      COMMON /FAURE/S,QS,COEF,NEXTN,TESTN,HISUM,RQS
      SAVE/FAURE/
C
C       FIND QUASI(1) USING FAURE (SECTION 3.3)
C
C       NEXTN HAS A REPRESENTATION IN BASE
C       QS OF THE FORM: SUM OVER J FROM ZERO
C       TO HISUM OF YTEMP(J)*(QS**J)
C
C       WE NOW COMPUTE THE YTEMP(J)'S.
C
      KTEMP = TESTN
      LTEMP = NEXTN
      DO 10 I = HISUM,0,-1
         KTEMP = KTEMP/QS
         MTEMP = MOD(LTEMP,KTEMP)
         YTEMP(I) = (LTEMP-MTEMP)/KTEMP
         LTEMP = MTEMP
   10 CONTINUE
C
C       QUASI(K) HAS THE FORM SUM OVER J
C       FROM ZERO TO HISUM OF
C       YTEMP(J)*(QS**(-(J+1)))
C
C       READY TO COMPUTE QUASI(1)
C       USING NESTED MULTIPLICATION
C
      R = YTEMP(HISUM)
      DO 20 I = HISUM - 1,0,-1
         R = YTEMP(I) + RQS*R
   20 CONTINUE
      QUASI(1) = R*RQS
C
C       FIND THE OTHER S-1 COMPONENTS
C       OF QUASI USING "FAURE" (SECTIONS
C       3.2 AND 3.3)
C
      DO 50 K = 2,S
         QUASI(K) = 0.0
         R = RQS
         DO 40 J = 0,HISUM
            ZTEMP = 0
            DO 30 I = J,HISUM
               ZTEMP = ZTEMP + COEF(I,J)*YTEMP(I)
C
C       NO APPARENT ALTERNATIVE
C       ONE-DIMENSIONAL COEFFICIENT ARRAY
C       EXCEPT VIA SUBSCRIPT ADDRESS
C       COMPUTATIONS AND EQUIVALENCING
C
   30       CONTINUE
C
C       NEW YTEMP(J) IS THE SUM
C       OVER I FROM J TO HISUM
C       OF (OLD YTEMP(I)*BINOM(I,J))
C       MOD QS
C
            YTEMP(J) = MOD(ZTEMP,QS)
            QUASI(K) = QUASI(K) + YTEMP(J)*R
            R = R*RQS
   40    CONTINUE
   50 CONTINUE
C
C       UPDATE NEXTN AND, IF NEEDED, TESTN AND
C       HISUM
C
      NEXTN = NEXTN + 1
      IF (NEXTN.EQ.TESTN) THEN
          TESTN = TESTN*QS
          HISUM = HISUM + 1
C
C       SINCE FLAG(2) IS TRUE,
C       HISUM STAYS UNDER 20
C
      END IF
C
      return
      end
      subroutine gohalt ( quasi )

c*********************************************************************72
c
cc GOHALT generates a new quasirandom Halton vector with each call.
C
C       IT ADAPTS KEY IDEAS FROM HALTON AND
C       SMITH, COMM. ACM. 7 (1964), 701-702.
C
C       THE USER MUST CALL "INHALT" BEFORE
C       CALLING "GOHALT". AFTER CALLING "INHALT",
C       TEST FLAG(1) AND FLAG(2); IF EITHER
C       IS FALSE, DO NOT CALL "GOHALT".
C       READ THE COMMENTS AT THE BEGINNING OF
C       "INHALT" AND THEN THOSE BELOW.
C
C       INPUTS:
C        FROM USER'S CALLING PROGRAM:
C        QUASI - THE LAST QUASIRANDOM VECTOR
C        GENERATED
C
C        FROM LABELLED COMMON "HALTON":
C        E           TOLERANCE
C        PRIME       PRIME(K) IS NOW THE RECIPROCAL OF
C                    THE K-TH PRIME NUMBER STARTING FROM TWO
C        S           DIMENSION
C
C       OUTPUT:
C        QUASI - A NEW QUASIRANDOM VECTOR
C
      DOUBLE PRECISION QUASI(40),PRIME(40),E,T,F,G,H

      INTEGER S,I
      COMMON /HALTON/E,PRIME,S
      SAVE/HALTON/
      EXTERNAL BDHALT
C
C  GENERATE QUASI ONE COMPONENT AT A TIME
C  USING RADIX PRIME(K) FOR COMPONENT K
C
      DO I = 1, S

         T = 1.0D+00 / PRIME(I)
         F = 1.0D+00 - QUASI(I)
         G = 1.0D+00
         H = T

   10    continue

         IF ( F - H .LT. E ) THEN
C
C  THIS CHECKS WHETHER Q+H>1-E
C
             G = H
             H = H * T
C
C  IF THIS IS THE (K-1)-ST TIME THIS
C  STATEMENT IS REACHED, CHECK WHETHER
C  QUASI(I)+R**(-K) > 1-E
C
             GO TO 10
C
C  THE BLOCK FROM STATEMENT 50 TO HERE IS A "LOOP WHILE"
C
         END IF
C
C  FOR THE APPROPRIATE I (DEPENDING ON HOW
C  MANY TIMES THE LOOP WHILE BLOCK ABOVE
C  IS EXECUTED), NOW ADD H**(I+1)+H**(I)-1
C  TO THE OLD QUASI(I) TO GET THE NEXT QUASI(I)
C
         QUASI(I) = G + H - F

      end do

      return
      end
      subroutine gosobl ( quasi )

c*********************************************************************72
C
cc GOSOBL generates a new quasirandom SOBOL vector with each call.
C
C     IT ADAPTS THE IDEAS OF ANTONOV AND SALEEV,
C     USSR COMPUT. MATHS. MATH. PHYS. 19 (1980),
C     252 - 256
C
C     THE USER MUST CALL "INSOBL" BEFORE CALLING
C     "GOSOBL".  AFTER CALLING "INSOBL", TEST
C     FLAG(1) AND FLAG(2);  IF EITHER IS FALSE,
C     DO NOT CALL "GOSOBL".  "GOSOBL" CHECKS
C     THAT THE USER DOES NOT MAKE MORE CALLS
C     THAN HE SAID HE WOULD : SEE THE COMMENTS
C     TO "INSOBL".
C
C     INPUTS:
C       FROM USER'S CALLING PROGRAM:
C         NONE
C
C       FROM LABELLED COMMON /SOBOL/:
C         V        TABLE OF DIRECTION NUMBERS
C         S        DIMENSION
C         MAXCOL   LAST COLUMN OF V TO BE USED
C         COUNT    SEQUENCE NUMBER OF THIS CALL
C         LASTQ    NUMERATORS FOR LAST VECTOR GENERATED
C         RECIPD   (1/DENOMINATOR) FOR THESE NUMERATORS
C
      INTEGER V(40,30),S,MAXCOL,COUNT,LASTQ(40)
      integer poly(40)
      INTEGER I,I2,L,EXOR
      REAL RECIPD,QUASI(40)

      COMMON /SOBOL/ V, S, MAXCOL, COUNT, POLY, LASTQ, RECIPD

      SAVE /SOBOL/
C
C  FIND THE POSITION OF THE RIGHT-HAND ZERO IN COUNT
C
      L = 0
      I = COUNT

   10 continue

      L = L + 1
      I2 = I / 2

      IF ( I .NE. 2 * I2 ) THEN
         I = I2
         GO TO 10
      END IF
C
C  CHECK THAT THE USER IS NOT CHEATING !
C
      IF ( MAXCOL .lt. l ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GOSOBL - Fatal error!'
        write ( *, '(a)' ) '  TOO MANY CALLS ON GOSOBL'
        stop
      end if
C
C  CALCULATE THE NEW COMPONENTS OF QUASI,
C  FIRST THE NUMERATORS, THEN NORMALIZED
C
      DO I = 1, S

c        LASTQ(I) = EXOR ( LASTQ(I), V(I,L) )
         LASTQ(I) = ieor ( LASTQ(I), V(I,L) )
C
C     IF A FULL-WORD EXCLUSIVE-OR, SAY .XOR., IS AVAILABLE,
C     THEN REPLACE THE PRECEDING STATEMENT BY
C
C       LASTQ(I) = LASTQ(I) .XOR. V(I,L)
C
C     TO GET A FASTER (EXTENDED FORTRAN) PROGRAM
C
         QUASI(I) = LASTQ(I) * RECIPD

      end do

      COUNT = COUNT + 1

      return
      end
      subroutine infaur ( flag, dimen, atmost, qss )

c*********************************************************************72
c
cc INFAUR initializes a Faure calculation.
c
C       THIS SUBROUTINE FIRST CHECKS WHETHER
C       THE USER-SUPPLIED DIMENSION "DIMEN" OF THE
C       QUASIRANDOM VECTORS IS ACCEPTABLE
C       (STRICTLY BETWEEN 1 AND 41) : IF SO,
C       FLAG(1)=.TRUE.
C
C       THEN IT CALCULATES AN UPPER SUMMATION
C       LIMIT "HISUM" BASED ON "DIMEN" AND THE
C       USER-SUPPLIED NUMBER "ATMOST" OF QUASIRANDOM
C       VECTORS REQUIRED. FLAG(2)=.TRUE. IF
C       ATMOST IS OK.
C
C       IF FLAG(1) AND FLAG(2) ARE TRUE,
C       "INFAUR" NEXT PRODUCES THE OTHER
C       OUTPUTS LISTED BELOW PASSED TO
C       SUBROUTINE GOFAUR VIA LABELLED
C       COMMON "FAURE". THESE OUTPUTS ARE
C       IRRELEVANT TO THE USER.
C
C       FIRST CALL INFAUR. IF FLAG(1) AND
C       FLAG(2) ARE TRUE, EACH (SUBSEQUENT)
C       CALL TO GOFAUR GENERATES A NEW
C       QUASIRANDOM VECTOR.
C
C       INPUTS : DIMEN, ATMOST
C
C       OUTPUTS
C          TO USERS CALLING PROGRAM:
C             FLAG
C             QSS   : SAME AS QS - SEE BELOW
C
C          TO GOFAUR:
C             S      :DIMENSION
C             QS     :SMALLEST PRIME >=S
C             COEF   :TABLE OF BINOMIAL
C                     COEFFICIENTS NEEDED
C                     BY GOFAUR.
C             NEXTN  :THE NUMBER OF THE
C                     NEXT QUASIRANDOM
C                     VECTOR,INITIALIZED
C                     TO TESTN-1 HERE.
C             TESTN  :INITIALIZED TO QS**4
C             HISUM  :AFTER BEING USED TO
C                     PRODUCE COEF, INITIALIZED
C                     TO 3 FOR GOFAUR.
C             RQS    :1.0/QS.
C
      LOGICAL FLAG(2)
C
      INTEGER S,ATMOST,QS,COEF(0:19,0:19),NEXTN,TESTN,HISUM,I,J,
     .        PRIMES(40),DIMEN
C
      REAL RQS
C
      COMMON /FAURE/S,QS,COEF,NEXTN,TESTN,HISUM,RQS
      SAVE/FAURE/
C
      DATA (PRIMES(I),I=1,40)/1,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,
     .     17,19,19,23,23,23,23,29,29,29,29,29,29,31,31,37,37,37,37,37,
     .     37,41,41,41/
C
C       CHECK S
C
      S = DIMEN
      FLAG(1) = S .GT. 1 .AND. S .LT. 41
      IF ( .NOT. FLAG(1)) RETURN
C
      QS = PRIMES(S)
      TESTN = QS**4
C
C         COMPUTE LOG(ATMOST+TESTN) IN BASE QS
C         USING A RATIO OF NATURAL LOGS TO GET
C         AN UPPER BOUND ON (THE NUMBER OF
C         DIGITS IN THE BASE QS REPRESENTATION
C         OF ATMOST+TESTN) MINUS ONE.
C
      HISUM = NINT(LOG(REAL(ATMOST+TESTN))/LOG(REAL(QS)))
      FLAG(2) = HISUM .LT. 20
      IF ( .NOT. FLAG(2)) RETURN
C
C        NOW FIND BINOMIAL COEFFICIENTS MOD QS
C        IN A LOWER-TRIANGULAR MATRIX "COEF"
C        USING RECURSION BINOM(I,J)=BINOM(I-1,J)
C        +BINOM(I-1,J-1) AND A=B+C IMPLIES MOD(A,D)=
C        MOD(MOD(B,D)+MOD(C,D),D)
C
      COEF(0,0) = 1
      DO 10 J = 1,HISUM
         COEF(J,0) = 1
         COEF(J,J) = 1
   10 CONTINUE
      DO 30 J = 1,HISUM
         DO 20 I = J + 1,HISUM
            COEF(I,J) = MOD(COEF(I-1,J)+COEF(I-1,J-1),QS)
   20    CONTINUE
   30 CONTINUE
C
C        CALCULATING THESE COEFFICIENTS
C        MOD QS AVOIDS POSSIBLE OVERFLOW
C        PROBLEMS WITH RAW BINOMIAL COEFFICIENTS
C
C        NOW COMPLETE INITIALIZATION
C        AS DESCRIBED IN SECTION 2.
C        NEXTN HAS 4 DIGITS IN BASE
C        QS, SO HISUM EQUALS 3.
C
      NEXTN = TESTN - 1
      HISUM = 3
      RQS = 1.0/REAL(QS)
C
      return
      end
      subroutine inhalt ( flag, dimen, atmost, tiny, quasi )

c*********************************************************************72
c
cc INHALT initializes a Halton calculation.
c
C       THIS SUBROUTINE FIRST CHECKS WHETHER
C       THE USER-SUPPLIED DIMENSION "DIMEN" OF
C       THE QUASIRANDOM VECTORS IS ACCEPTABLE
C       (STRICTLY BETWEEN 0 AND 41): IF SO,
C       FLAG(1)=.TRUE.
C
C       THE USER SUPPLIES
C         ATMOST : AN UPPER BOUND ON
C                  THE NUMBERS OF VECTORS
C                  REQUIRED
C         TINY   : 2**(-B), WHERE B IS THE NUMBER
C                  OF BITS IN THE MANTISSAS OF
C                  DOUBLE-PRECISION CONSTANTS, OR
C                  SLIGHTLY HIGHER
C
C       IF FLAG(1) IS TRUE, THEN INHALT
C       CALCULATES A TOLERANCE PARAMETER
C       "E" TO MAKE THE PROGRAM WORK
C       CORRECTLY IN FINITE PRECISION
C       ARITHMETIC AND A PARAMETER "DELTA"
C       TO CHECK THAT "E" WORKS. IF ALL IS
C       WELL, FLAG(2) GETS SET TO TRUE.
C       IF FLAG(2) IS FALSE, THEN ATMOST
C       IS TOO BIG RELATIVE TO TINY.
C
C       IF FLAG(1) AND FLAG(2) ARE TRUE,
C       INHALT COMPUTES THE FIRST VECTOR
C       "QUASI". ALL SUBSEQUENT VECTORS
C       COME FROM GOHALT.
C
C       INPUTS: DIMEN,ATMOST,TINY,PRIME
C
C       GET PRIME FROM LABELLED COMMON
C       "HALTON" INITIALIZED WITH
C       BLOCK DATA "BDHALT"
C
C       OUTPUTS:
C         TO USER'S CALLING PROGRAM: FLAG,QUASI
C
C         TO GOHALT VIA LABELLED COMMON
C         "HALTON":
C            E     : TOLERANCE
C            PRIME : PRIME(K) IS THE RECIPROCAL
C                    OF THE K-TH PRIME NUMBER
C                    STARTING FROM TWO
C            S     : DIMENSION
C
      LOGICAL FLAG(2)
C
      INTEGER S,ATMOST,I,DIMEN
C
      DOUBLE PRECISION QUASI(40),PRIME(40),E,TINY,DELTA
C
      COMMON /HALTON/E,PRIME,S
      SAVE/HALTON/
      EXTERNAL BDHALT
C
C  CHECK DIMEN
C
      S = DIMEN

      FLAG(1) = S .GT. 1 .AND. S .LT. 41

      IF ( .NOT. FLAG(1)) then
        RETURN
      end if
C
C  COMPUTE AND CHECK TOLERANCE
C
      E = 0.9D+00 
     &  * ( 1.0D+00 / ( ATMOST * PRIME ( S ) ) - 10.0D+00 * TINY )

      DELTA = 100 * TINY * DBLE ( ATMOST + 1 ) 
     &  * DLOG10 ( DBLE ( ATMOST ) )

      FLAG(2) = DELTA .LE. 0.09D+00 * ( E - 10.0D+00 * TINY )

      IF ( .NOT. FLAG(2) ) then
        RETURN
      end if
C
C  COMPUTE FIRST VECTOR.
C
      DO I = 1, S
        QUASI(I) = 1.0D+00 / PRIME(I)
      end do

      return
      end
      subroutine insobl ( flag, dimen, atmost, taus )

c*********************************************************************72
c
cc INSOBL initializes a Sobol calculation.
c
C     FIRST CHECK WHETHER THE USER-SUPPLIED
C     DIMENSION "DIMEN" OF THE QUASI-RANDOM
C     VECTORS IS STRICTLY BETWEEN 0 AND 41.
C     IF SO, FLAG(1) = .TRUE.
C
C     NEXT CHECK "ATMOST", AN UPPER BOUND ON THE NUMBER
C     OF CALLS THE USER INTENDS TO MAKE ON "GOSOBL".  IF
C     THIS IS LESS THAN 2**30, THEN FLAG(2) = .TRUE.
C     (WE ASSUME WE ARE WORKING ON A COMPUTER WITH
C     WORD LENGTH AT LEAST 31 BITS EXCLUDING SIGN.)
C     THE NUMBER OF COLUMNS OF THE ARRAY V WHICH
C     ARE INITIALIZED IS
C          MAXCOL = NUMBER OF BITS IN ATMOST.
C     IN "GOSOBL" WE CHECK THAT THIS IS NOT EXCEEDED.
C
C     THE LEADING ELEMENTS OF EACH ROW OF V ARE
C     INITIALIZED BY "BDSOBL".
C     EACH ROW CORRESPONDS TO A PRIMITIVE POLYNOMIAL
C     (AGAIN, SEE "BDSOBOL").  IF THE POLYNOMIAL HAS
C     DEGREE M, ELEMENTS AFTER THE FIRST M ARE CALCULATED.
C
C     THE NUMBERS IN V ARE ACTUALLY BINARY FRACTIONS.
C     "RECIPD" HOLDS 1/(THE COMMON DENOMINATOR OF ALL
C     OF THEM).
C
C     "INSOBL" IMPLICITLY COMPUTES THE FIRST ALL-ZERO
C     VECTOR, BUT DOES NOT RETURN IT TO THE CALLING
C     PROGRAM. SUBSEQUENT VECTORS COME FROM "GOSOBL".
C     AFTER INITIALIZATION IS COMPLETE, WE RE-USE THE SPACE
C     OCCUPIED BY "POLY" TO HOLD "LASTQ", THE
C     NUMERATORS OF THE LAST VECTOR GENERATED.
C
C     "TAUS" IS FOR DETERMINING "FAVORABLE" VALUES. AS
C     DISCUSSED IN BRATLEY/FOX, THESE HAVE THE FORM
C     N = 2**K WHERE K .GE. (TAUS+S-1) FOR INTEGRATION
C     AND K .GT. TAUS FOR GLOBAL OPTIMIZATION.
C
C     INPUTS :
C       FROM USER'S PROGRAM : DIMEN, ATMOST
C       FROM BLOCK DATA "BDSOBL" : POLY, PART OF V
C
C     OUTPUTS :
C       TO USER'S PROGRAM : FLAG, TAUS
C       TO "GOSOBL" VIA /SOBOL/ :
C         V, COUNT, LASTQ, MAXCOL, S
C
      integer atmost
      integer count
      integer dimen
      integer exor
      logical flag(2)
      integer i
      logical includ(8)
      integer j
      integer j2
      integer k
      integer l
      integer lastq(40)
      integer m
      integer maxcol
      integer newv
      integer poly(40)
      real recipd
      integer s
      integer tau(13)
      integer taus
      integer v(40,30)

      COMMON /SOBOL/ V, S, MAXCOL, COUNT, POLY, LASTQ, RECIPD

      SAVE /SOBOL/
      save tau

      DATA TAU / 0,0,1,3,5,8,11,15,19,23,27,31,35 /
C
C  CHECK PARAMETERS
C
      S = DIMEN

      FLAG(1) = (S.GE.1 .AND. S.LE.40)

      FLAG(2) = ( ATMOST .LT. 2**30 )

      IF ( .NOT. (FLAG(1).AND.FLAG(2))) RETURN

      IF ( S .LE. 13 ) THEN
        TAUS = TAU(S)
      ELSE
        TAUS = -1
      END IF
C
C  FIND NUMBER OF BITS IN ATMOST
C
      I = ATMOST
      MAXCOL = 0

   10 continue

      MAXCOL = MAXCOL + 1
      I = I / 2

      IF ( 0 .lt. I ) then
        GO TO 10
      end if
C
C  INITIALIZE ROW 1 of V.
C
      DO I = 1, MAXCOL
        V(1,I) = 1
      end do
C
C  INITIALIZE REMAINING ROWS of V.
C
      DO I = 2, S
C
C  THE BIT PATTERN OF POLYNOMIAL I GIVES ITS FORM
C  (SEE COMMENTS TO "BDSOBL")
C  FIND DEGREE OF POLYNOMIAL I FROM BINARY ENCODING
C
         J = POLY(I)
         M = 0

   30    continue

         J = J / 2

         IF ( 0 .lt. J ) THEN
             M = M + 1
             GO TO 30
         END IF
C
C  WE EXPAND THIS BIT PATTERN TO SEPARATE COMPONENTS
C  OF THE LOGICAL ARRAY INCLUD.
C
         J = POLY(I)

         DO K = M,1,-1
            J2 = J/2
            INCLUD(K) = (J.NE.2*J2)
            J = J2
         end do
C
C  CALCULATE REMAINING ELEMENTS OF ROW I AS EXPLAINED
C  IN BRATLEY AND FOX, SECTION 2
C
         DO J = M + 1,MAXCOL
            NEWV = V(I,J-M)
            L = 1
            DO K = 1,M
               L = 2*L
               IF (INCLUD(K)) then
                 NEWV = EXOR(NEWV,L*V(I,J-K))
               end if
C
C     IF A FULL-WORD EXCLUSIVE-OR, SAY .XOR., IS AVAILABLE,
C     THEN REPLACE THE PRECEDING STATEMENT BY
C
C           IF (INCLUD(K)) NEWV = NEWV .XOR. L * V(I, J-K)
C
C     TO GET A FASTER (EXTENDED FORTRAN) PROGRAM
C
            end do
            V(I,J) = NEWV
         end do

      end do
C
C  MULTIPLY COLUMNS OF V BY APPROPRIATE POWER OF 2
C
      L = 1
      DO J = MAXCOL - 1, 1, -1
        L = 2 * L
        DO I = 1, S
          V(I,J) = V(I,J) * L
        end do
      end do
C
C  RECIPD IS 1/(COMMON DENOMINATOR OF THE ELEMENTS IN V)
C
      RECIPD = 1.0E+00 / real ( 2 * L )
C
C  SET UP FIRST VECTOR AND VALUES FOR "GOSOBL"
C
      COUNT = 0
      DO I = 1, S
        LASTQ(I) = 0
      end do

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
      function unif ( ix )

c*********************************************************************72
c
cc UNIF is a portable pseudorandom number generator.
c
C       PORTABLE PSEUDORANDOM NUMBER
C       GENERATOR IMPLEMENTING THE RECURSION
C
C       IX=16807*IX MOD(2**31-1)
C       UNIF=IX/(2**31-1)
C
C       USING ONLY 32 BITS INCLUDING
C       SIGN
C
C       INPUT:
C        IX =INTEGER STRICTLY BETWEEN
C           0 AND 2** 31 -1
C
C       OUTPUTS:
C        IX=NEW PSEUDORANDOM INTEGER
C           STRICTLY BETWEEN 0 AND
C           2**31-1
C        UNIF=UNIFORM VARIATE (FRACTION)
C             STRICTLY BETWEEN 0 AND 1
C
C      FOR JUSTIFICATION, SEE P. BRATLEY,
C      B.L. FOX, AND L.E. SCHRAGE (1983)
C      "A GUIDE TO SIMULATION"
C      SPRINGER-VERLAG, PAGES 201-202
C
      INTEGER K1,IX
      real unif

      K1 = IX / 127773
      IX = 16807 * (IX-K1*127773) - K1*2836
      IF (IX.LT.0) IX = IX + 2147483647
      UNIF = IX * 4.656612875E-10

      RETURN
      END
