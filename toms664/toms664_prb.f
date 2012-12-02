      program main

c*********************************************************************72
C
C     TESTDRIVER FOR THE SUBROUTINE GBSOL               JANUARY 15, 1988
C
C     THIS PROGRAM IS A 1977 AMERICAN NATIONAL STANDARD FORTRAN IMPLE-
C     MENTATION OF A TESTDRIVER FOR THE SUBROUTINE GBSOL.
C
C     THERE ARE SIX TESTS BASED ON THE FOLLOWING BAND MATRIX:
C
C
C      20.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0.  0.  0.  0. ...
C
C      -1.  0.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0.  0.  0. ...
C
C       2. -1. 20.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0.  0. ...
C
C      -3.  2. -1.  0.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0. ...
C
C       4. -3.  2. -1. 20.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10. ...
C
C      -5.  4. -3.  2. -1.  0.  1. -2.  3. -4.  5. -6.  7. -8.  9. ...
C
C       6. -5.  4. -3.  2. -1. 20.  1. -2.  3. -4.  5. -6.  7. -8. ...
C
C      -7.  6. -5.  4. -3.  2. -1.  0.  1. -2.  3. -4.  5. -6.  7. ...
C
C       8. -7.  6. -5.  4. -3.  2. -1. 20.  1. -2.  3. -4.  5. -6. ...
C
C      -9.  8. -7.  6. -5.  4. -3.  2. -1.  0.  1. -2.  3. -4.  5. ...
C
C      10. -9.  8. -7.  6. -5.  4. -3.  2. -1. 20.  1. -2.  3. -4. ...
C
C       0. 10. -9.  8. -7.  6. -5.  4. -3.  2. -1.  0.  1. -2.  3. ...
C
C       0.  0. 10. -9.  8. -7.  6. -5.  4. -3.  2. -1. 20.  1. -2. ...
C
C       0.  0.  0. 10. -9.  8. -7.  6. -5.  4. -3.  2. -1.  0.  1. ...
C
C       0.  0.  0.  0. 10. -9.  8. -7.  6. -5.  4. -3.  2. -1. 20. ...
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C
C
C     THE ROWS OF THIS MATRIX ARE CALCULATED BY THE SUBROUTINE ROW WHICH
C     IS SUPPLIED TOGETHER WITH THE TESTDRIVER.
C
C     IN ALL TESTS, THE RIGHT HAND SIDE IS CHOSEN SO THAT ALL COMPONENTS
C     OF THE SOLUTION VECTOR HAVE THE VALUE ONE.
C
C     IN THE FIRST FOUR TESTS A SYSTEM WITH ONE THOUSAND EQUATIONS IS
C     SOLVED.  TO BEGIN WITH, THE SECOND DIMENSION NA OF A IS CHOSEN AS
C     1000 SO THAT NO I/O OPERATIONS ARE PERFORMED. THEN, NA IS SET TO
C     THE VALUE OF 511 SO THAT THE MATRIX IS DIVIDED INTO TWO BLOCKS AND
C     THE FORTRAN 77 VERSIONS OF THE SUBROUTINES GOPEN, GCLOSE, AND
C     GRWRAN ARE TESTED (THE USER SHOULD REWRITE THESE SUBROUTINES BY
C     APPLYING ASSEMBLER UTILITY PROGRAMS OF HIS INSTALLATION. FORTRAN
C     SUPPORTED I/O IS TOO INEFFICIENT AND SHOULD BE USED ONLY FOR
C     TESTING). IN THE THIRD TEST, THE VALUE OF NA IS TOO SMALL, AND
C     GBSOL RETURNS TO THE CALLING PROGRAM WITH THE ERROR CODE "JERR=1."
C     IN THE FOURTH TEST THE PARAMETER EPS IS CHOSEN TO BE SLIGHTLY
C     LARGER THAN THE SMALLEST PIVOT SO THAT THE PIVOT IS REPLACED BY
C     EPS.
C
C     FOR THE LAST TWO TESTS, THE FIRST LINE OF THE MATRIX IS CHANGED
C     IN ORDER TO PRODUCE A NEARLY SINGULAR MATRIX WITH A SMALL PIVOT:
C
C
C       0.  0.  0.  0.  0.  0.  0.  0.  0.  0. DEL  0.  0.  0.  0. ...
C
C      -1.  0.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0.  0.  0. ...
C
C       2. -1. 20.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0.  0. ...
C
C      -3.  2. -1.  0.  1. -2.  3. -4.  5. -6.  7. -8.  9.-10.  0. ...
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C       .   .   .   .   .   .   .   .   .   .   .   .   .   .   .
C
C
C     DEL IS SET TO THE VALUE OF 1.D-12.
C
C     IN THE FIFTH TEST, THE DIMENSION N IS 750, NA IS SET TO 241, AND
C     THE SYSTEM IS SOLVED. FINALLY, THE PARAMETER EPS IS AGAIN CHOSEN
C     TO BE SLIGHTLY LARGER THAN THE SMALLEST PIVOT SO THAT THE PIVOT
C     IS REPLACED BY EPS.
C
C     AUTHOR: GEZA SCHRAUF, CAL TECH
C
C***********************************************************************
      PARAMETER (NP=1000,NAP=1000)

      DOUBLE PRECISION A(21,NAP),B(NP,1)
      DOUBLE PRECISION ALOGDT,DEL,EPS,ERROR,PIVMAX,PIVMIN,SIGNDT
      INTEGER IFLAG,INFO,IPIV,JERR,L,M,MD,N,NA,NFSUA,NFSUAB,NOEOA,
     .        NOEOB,NREAL8,NRHS,N1

      COMMON /TOROW/ DEL,IFLAG,N

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS664_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS664 library.'

      IFLAG = 1
      M = 21
      MD = (M+1)/2
      N = NP
      NA = NAP
      NRHS = 1
      N1 = N - 1
C-----------------------------------------------------------------------
C               FIRST TEST:  SOLVE THE SYSTEM WITHOUT DIRECT ACCESS I/O.
C-----------------------------------------------------------------------
      WRITE(6,9001)
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
      NOEOA = M*NA
      NFSUA = 2*NOEOA
      NREAL8 = NOEOA + N
      NFSUAB = 2*NREAL8
C
      WRITE (6,9021) NFSUAB,NFSUA
C
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
C
      DO 10 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
   10 CONTINUE
C
      B(1,1) = 15.D0
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 1.D-12
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9041) ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,INFO,JERR
C
C                                      CALCULATE THE ERROR IN L2 - NORM.
      ERROR = 0.D0
C
      DO 20 L = 1,N
         ERROR = ERROR + (B(L,1)-1.D0)**2
   20 CONTINUE
C
      ERROR = DSQRT(ERROR/N)
C
      WRITE (6,9051) ERROR
C                                 CALCULATE THE ERROR IN MAXIMUM - NORM.
      ERROR = 0.D0
C
      DO 30 L = 1,N
         TEMP = DABS(B(L,1)-1.D0)
         IF (TEMP.GT.ERROR) ERROR = TEMP
   30 CONTINUE
C
      WRITE (6,9052) ERROR
C-----------------------------------------------------------------------
C                SECOND TEST:   SOLVE THE SYSTEM WITH DIRECT ACCESS I/O.
C-----------------------------------------------------------------------
      WRITE(6,9002)
C
      NA = 511
C
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
      NOEOA = M*NA
      NFSUA = 2*NOEOA
      NREAL8 = NOEOA + N
      NFSUAB = 2*NREAL8
C
      WRITE (6,9021) NFSUAB,NFSUA
C
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
C
      DO 40 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
   40 CONTINUE
C
      B(1,1) = 15.D0
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 1.D-12
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9041) ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,INFO,JERR
C
C                                      CALCULATE THE ERROR IN L2 - NORM.
      ERROR = 0.D0
C
      DO 50 L = 1,N
         ERROR = ERROR + (B(L,1)-1.D0)**2
   50 CONTINUE
C
      ERROR = DSQRT(ERROR/N)
C
      WRITE (6,9051) ERROR
C                                 CALCULATE THE ERROR IN MAXIMUM - NORM.
      ERROR = 0.D0
C
      DO 60 L = 1,N
         TEMP = DABS(B(L,1)-1.D0)
         IF (TEMP.GT.ERROR) ERROR = TEMP
   60 CONTINUE
C
      WRITE (6,9052) ERROR
C-----------------------------------------------------------------------
C                                THIRD TEST:   VALUE OF NA IS TOO SMALL.
C-----------------------------------------------------------------------
      WRITE(6,9003)
C
      NA = 7
C
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
C
      DO 70 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
   70 CONTINUE
C
      B(1,1) = 15.D0
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 1.D-12
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9043) JERR
C
C-----------------------------------------------------------------------
C            FOURTH TEST:   SOLVE A SYSTEM WITH REPLACED SMALLEST PIVOT.
C-----------------------------------------------------------------------
      WRITE(6,9004)
C
      NA = 511
C
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
      NOEOA = M*NA
      NFSUA = 2*NOEOA
      NREAL8 = NOEOA + N
      NFSUAB = 2*NREAL8
C
      WRITE (6,9021) NFSUAB,NFSUA
C
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
C
      DO 80 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
   80 CONTINUE
C
      B(1,1) = 15.D0
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 0.9667910D+1
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9041) ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,INFO,JERR
C
C                                      CALCULATE THE ERROR IN L2 - NORM.
      ERROR = 0.D0
C
      DO 90 L = 1,N
         ERROR = ERROR + (B(L,1)-1.D0)**2
   90 CONTINUE
C
      ERROR = DSQRT(ERROR/N)
C
      WRITE (6,9051) ERROR
C                                 CALCULATE THE ERROR IN MAXIMUM - NORM.
      ERROR = 0.D0
C
      DO 100 L = 1,N
         TEMP = DABS(B(L,1)-1.D0)
         IF (TEMP.GT.ERROR) ERROR = TEMP
  100 CONTINUE
C
      WRITE (6,9052) ERROR
C-----------------------------------------------------------------------
C                       FIFTH TEST:   SOLVE A SYSTEM WITH A SMALL PIVOT.
C-----------------------------------------------------------------------
      WRITE(6,9005)
C
      IFLAG = 2
C
      N = 750
      NA = 241
C
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
      NOEOA = M*NA
      NFSUA = 2*NOEOA
      NREAL8 = NOEOA + N
      NFSUAB = 2*NREAL8
C
      WRITE (6,9021) NFSUAB,NFSUA
C
      DEL = 1.D-12
      WRITE (6,9031) DEL
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
      DO 110 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
  110 CONTINUE
C
      B(1,1) = DEL
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 1.D-12
      WRITE (6,9032) EPS
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9041) ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,INFO,JERR
C
C                                      CALCULATE THE ERROR IN L2 - NORM.
      ERROR = 0.D0
C
      DO 120 L = 1,N
         ERROR = ERROR + (B(L,1)-1.D0)**2
  120 CONTINUE
C
      ERROR = DSQRT(ERROR/N)
C
      WRITE (6,9051) ERROR
C
C                                 CALCULATE THE ERROR IN MAXIMUM - NORM.
      ERROR = 0.D0
C
      DO 130 L = 1,N
         TEMP = DABS(B(L,1)-1.D0)
         IF (TEMP.GT.ERROR) ERROR = TEMP
  130 CONTINUE
C
      WRITE (6,9052) ERROR
C-----------------------------------------------------------------------
C             SIXTH TEST:   SOLVE A SYSTEM WITH REPLACED SMALLEST PIVOT.
C-----------------------------------------------------------------------
      WRITE(6,9006)
C
      IFLAG = 2
C
      N = 750
      NA = 241
C
      NOEOB = M*N
      WRITE (6,9011) M,NA,N,NOEOB
C
      NOEOA = M*NA
      NFSUA = 2*NOEOA
      NREAL8 = NOEOA + N
      NFSUAB = 2*NREAL8
C
      WRITE (6,9021) NFSUAB,NFSUA
C
      DEL = 1.D-12
      WRITE (6,9031) DEL
C                                  ASSIGN VALUES TO THE RIGHT HAND SIDE.
      DO 140 L = 1,N,2
         B(L,1) = 20.D0
         B(L+1,1) = 0.D0
  140 CONTINUE
C
      B(1,1) = DEL
      B(2,1) = -6.D0
      B(3,1) = 16.D0
      B(4,1) = -7.D0
      B(5,1) = 17.D0
      B(6,1) = -8.D0
      B(7,1) = 18.D0
      B(8,1) = -9.D0
      B(9,1) = 19.D0
      B(10,1) =-10.D0
C
      B(N-9,1) = 30.D0
      B(N-8,1) =  1.D0
      B(N-7,1) = 29.D0
      B(N-6,1) =  2.D0
      B(N-5,1) = 28.D0
      B(N-4,1) =  3.D0
      B(N-3,1) = 27.D0
      B(N-2,1) =  4.D0
      B(N-1,1) = 26.D0
      B(N,1) = 5.D0
C
      EPS = 0.22837D-9
      WRITE (6,9032) EPS
      INFO = 0
C
C                                            SOLVE THE SYSTEM   A*X = B.
C
      CALL GBSOL(A,M,NA,B,N,NRHS,ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,EPS,
     .           INFO,JERR)
C
      WRITE (6,9041) ALOGDT,SIGNDT,PIVMAX,PIVMIN,IPIV,INFO,JERR
C
C                                      CALCULATE THE ERROR IN L2 - NORM.
      ERROR = 0.D0
C
      DO 150 L = 1,N
         ERROR = ERROR + (B(L,1)-1.D0)**2
  150 CONTINUE
C
      ERROR = DSQRT(ERROR/N)
C
      WRITE (6,9051) ERROR
C                                 CALCULATE THE ERROR IN MAXIMUM - NORM.
      ERROR = 0.D0
C
      DO 160 L = 1,N
         TEMP = DABS(B(L,1)-1.D0)
         IF (TEMP.GT.ERROR) ERROR = TEMP
  160 CONTINUE
C
      WRITE (6,9052) ERROR
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS664:'
      write ( *, '(a)' ) '  Normal end of execution.'

      STOP
C
 9001 FORMAT(//' ',71('-')/' F I R S T   T E S T   W I T H O U T  ',
     .       ' D I R E C T - A C C E S S   I / O'/' ',71('-')//)
 9002 FORMAT(//' ',68('-')/' S E C O N D   T E S T   W I T H   ',
     .       ' D I R E C T - A C C E S S   I / O'/' ',68('-')//)
 9003 FORMAT(//' ',71('-')/' T H I R D   T E S T   W I T H   T O O  ',
     .       ' S M A L L   V A L U E   OF   NA'/' ',71('-')//)
 9004 FORMAT(//' ',79('-')/' F O U R T H   T E S T   W I T H  ',
     .       ' S M A L L E S T   P I V O T   R E P L A C E D'/
     .       ' ',79('-')//)
 9005 FORMAT(//' ',69('-')/' F I F T H   T E S T   W I T H  ',
     .       ' S M A L L   P I V O T   E L E M E N T'/' ',69('-')//)
 9006 FORMAT(//' ',77('-')/' S I X T H   T E S T   W I T H  ',
     .       ' S M A L L E S T   P I V O T   R E P L A C E D'/
     .       ' ',77('-')//)
 9011 FORMAT (//' FIRST DIMENSION OF A:                            ',
     .  ' M     =',I10/
     .  ' SECOND DIMENSION OF A:                           ',
     .  ' NA    =',I10//
     .  ' NUMBER OF EQUATIONS:                             ',
     .  ' N     =',I10/
     .  ' NUMBER OF ELEMENTS IN THE BAND OF THE MATRIX A:  ',
     .  ' M*N   =',I10///)
 9021 FORMAT (' MEMORY REQUIREMENTS:'//
     .  ' NUMBER OF FORTRAN STORAGE UNITS NEEDED     '/
     .  ' TO STORE THE MATRIX A AND THE RHS B:        ',
     .  ' 2*(M*NA+N) =',I10//
     .  ' NUMBER OF FORTRAN STORAGE UNITS NEEDED     '/
     .  ' TO STORE THE MATRIX A:                      ',
     .  '   2*(M*NA) =',I10//)
 9031 FORMAT (' DEL    = ',E14.7/)
 9032 FORMAT (' EPS    = ',E14.7)
 9041 FORMAT (//' ALOGDT = ',E14.7/' SIGNDT = ',F3.0//' PIVMAX = ',
     .  E14.7/' PIVMIN = ',E14.7/' IPIV   = ',I10//' INFO   = ',
     .  I10//' JERR   = ',I10)
 9043 FORMAT (//' JERR   = ',I10//
     .          ' RETURN TO CALLING PROGRAM WITH ERROR CODE',
     .          ' "JERR = 1,"'/' BECAUSE  NA = 7 ',
     .          ' IS LESS THAN  (M+1)/2 + 2 = 13.')
 9051 FORMAT (//' ERROR IN THE L2-NORM:   ',E15.7)
 9052 FORMAT (' AND IN THE MAXIMUM-NORM:',E15.7//)
      END
      SUBROUTINE ROW(AROW,L)

c*********************************************************************72

      DOUBLE PRECISION AROW(*),DEL,ATEMP
      INTEGER I,IFLAG,ISIGN,ITEMP,L,MD,MDMI,MDPI,N
C
      COMMON /TOROW/ DEL,IFLAG,N
C
      ISIGN = -1
      MD = 11
C
      ITEMP = L/2
      ATEMP = L - ITEMP*2
C
      AROW(MD) = 20.D0 * ATEMP
C
      DO 10 I = 1,10
         ISIGN = ISIGN* (-1)
         MDMI = MD - I
         MDPI = MD + I
C
         AROW(MDMI) = -ISIGN*I
         AROW(MDPI) = ISIGN*I
   10 CONTINUE
C
      IF (L.GT.10) GO TO 30
C
      IF (L.LE.10) AROW(1) = 0.D0
      IF (L.LE.9) AROW(2) = 0.D0
      IF (L.LE.8) AROW(3) = 0.D0
      IF (L.LE.7) AROW(4) = 0.D0
      IF (L.LE.6) AROW(5) = 0.D0
      IF (L.LE.5) AROW(6) = 0.D0
      IF (L.LE.4) AROW(7) = 0.D0
      IF (L.LE.3) AROW(8) = 0.D0
      IF (L.LE.2) AROW(9) = 0.D0
      IF (L.GT.1) GO TO 30
      AROW(10) = 0.D0
      IF (IFLAG.EQ.1) GO TO 30
      DO 20 I=11,20
         AROW(I) = 0.D0
   20 CONTINUE
      AROW(21) = DEL
C
   30 IF (L.LT.N-9) GO TO 40
C
      IF (L.GE.N) AROW(12) = 0.D0
      IF (L.GE.N-1) AROW(13) = 0.D0
      IF (L.GE.N-2) AROW(14) = 0.D0
      IF (L.GE.N-3) AROW(15) = 0.D0
      IF (L.GE.N-4) AROW(16) = 0.D0
      IF (L.GE.N-5) AROW(17) = 0.D0
      IF (L.GE.N-6) AROW(18) = 0.D0
      IF (L.GE.N-7) AROW(19) = 0.D0
      IF (L.GE.N-8) AROW(20) = 0.D0
      IF (L.GE.N-9) AROW(21) = 0.D0
   40 CONTINUE
C
      RETURN
      END
