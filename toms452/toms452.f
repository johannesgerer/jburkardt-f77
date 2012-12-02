      SUBROUTINE NXCBN(N, M, IC)

c*********************************************************************72
c
C  EXPLANATION OF THE PARAMETERS IN THE CALLING SEQUENCE:
C    N, THE TOTAL NUMBER OF OBJECTS.
C    M, THE NUMBER OF OBJECTS TO BE TAKEN FROM N.
C      IF M = 0, OR M>=N, EXIT WITH ARGUMENTS UNCHANGED.
C    IC, AN INTEGER ARRAY.  IC CONTAINS AN N-DIMNSIONAL
C      BINARY VECTOR WITH M ELEMENTS SET TO 1, REPRESENTING
C      THE N OBJECTS IN A COMBINATION.
C  THIS ALGORITHM IS PROGRAMMED IN ANSI STANDARD FORTRAN.
c
      INTEGER IC(N)
c
C  CHECK ENDING PATTERN OF VECTOR.
c
      IF (M.GE.N .OR. M.EQ.0 ) GO TO 140
      N1 = N - 1
      DO 10 J=1,N1
        NJ = N - J
        IF (IC(N).EQ.IC(NJ)) GO TO 10
        J1 = J
        GO TO 20
10    CONTINUE
20    IF (MOD(M,2).EQ.1) GO TO 90
c
C  FOR M EVEN.
c
      IF ( IC(N).EQ.1 ) GO TO 30
      K1 = N - J1
      K2 = K1 + 1
      GO TO 130
30    IF (MOD(J1,2).EQ.1) GO TO 40
      GO TO 120
c
C  SCAN FROM RIGHT TO LEFT.
c
40    JP = ( N - J1 ) - 1
      DO 50 I=1,JP
        I1 = JP + 2 - I
        IF(IC(I1).EQ.0) GO TO 50
        IF(IC(I1-1).EQ.1) GO TO 70
        GO TO 80
50    CONTINUE
60    K1 = 1
      K2 = ( N + 1 ) - M
      GO TO 130
70    K1 = I1 - 1
      K2 = N - J1
      GO TO 130
80    K1 = I1 - 1
      K2 = ( N + 1 ) - J1
      GO TO 130
c
C  FOR M ODD.
c
90    IF (IC(N).EQ.1) GO TO 110
      K2 = ( N - J1 ) - 1
      IF ( K2 .EQ. 0 ) GO TO 60
      IF ( IC(K2+1) .EQ. 1 .AND. IC(K2) .EQ. 1 ) GO TO 100
      K1 = K2 + 1
      GO TO 130
100   K1 = N
      GO TO 130
110   IF ( MOD(J1,2) .EQ. 1 ) GO TO 120
      GO TO 40
120   K1 = N - J1
      K2 = MIN0 ( ( K1 + 2 ), N )
c
C  COMPLEMENTING TWO BITS TO OBTAIN THE NEXT COMBINATION.
c
130   IC(K1) = 1 - IC(K1)
      IC(K2) = 1 - IC(K2)
140   RETURN
      END
