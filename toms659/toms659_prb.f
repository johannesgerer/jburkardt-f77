      PROGRAM TESTS

c*********************************************************************72
C
C       THIS PROGRAM TESTS ACCURACY OF
C       NUMERICAL INTEGRATION USING "GOSOBL"
C       AND INTEGRAND (2) OF DAVIS AND
C       RABINOWITZ, PAGE 406
C
C       IT USES A NONSTANDARD TIMING
C       ROUTINE "SECOND"
C
C       PARAMETER STATEMENT SPECIFIES INPUT
C       AND OUTPUT UNITS
C
C
C     .. Parameters ..
      INTEGER INPUT,OUTPUT
      PARAMETER (INPUT=5,OUTPUT=6)
      INTEGER MAXDIM
      PARAMETER (MAXDIM=1111)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION F,SUM
      INTEGER ATMOST,DIMEN,I,J,TAUS
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION QUASI(MAXDIM)
      LOGICAL FLAG(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL GOSOBL,INSOBL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
   10 READ (*,*) DIMEN,ATMOST
      IF (DIMEN.EQ.0) STOP ' RUN ENDS NORMALLY'
      WRITE (OUTPUT,FMT='(1H1)')
      WRITE (OUTPUT,FMT=*) 'TEST SOBOL'
      WRITE (OUTPUT,FMT=*) 'DIMENSION = ',DIMEN
      WRITE (OUTPUT,FMT=*) 'ATMOST = ',ATMOST
C
C     CALL SECOND(T1)
      CALL INSOBL(FLAG,DIMEN,ATMOST,TAUS)
      IF (.NOT.FLAG(1)) THEN
          WRITE (OUTPUT,FMT=*) 'DIMENSION = ',DIMEN
          WRITE (OUTPUT,FMT=*) 'DIMEN IS NOT OK'
          STOP

      END IF

      IF (.NOT.FLAG(2)) THEN
          WRITE (OUTPUT,FMT=*) 'ATMOST = ',ATMOST
          WRITE (OUTPUT,FMT=*) 'ATMOST IS NOT OK'
          STOP

      END IF
C     WRITE(OUTPUT,*) 'START TIME = ',T1
      WRITE (OUTPUT,FMT=*) 'I = ITERATION NUMBER'
      WRITE (OUTPUT,FMT=*) 'EI = ESTIMATE OF INTEGRAL'
      WRITE (OUTPUT,FMT='(1H )')
C
      SUM = 0.0
      DO 30 I = 1,ATMOST
          CALL GOSOBL(QUASI)
          F = 1.0
          DO 20 J = 1,DIMEN
              F = F*ABS(4.0*QUASI(J)-2.0)
   20     CONTINUE
          SUM = SUM + F
C         IF(MOD(I,500).EQ.0) THEN
C             WRITE(OUTPUT,*) 'I = ',I
C             WRITE(OUTPUT,*) 'EI = ',SUM/I
C             CALL SECOND(T2)
C             WRITE(OUTPUT,*) 'TIMEI = ',T2-T1
C             WRITE(OUTPUT,'(1H )')
C           ENDIF
   30 CONTINUE
      WRITE (OUTPUT,FMT=*) ' EI = ',SUM/ATMOST
C
      GO TO 10

      END
