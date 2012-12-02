      SUBROUTINE ROUTNG ( N, P, SN, EN, M, D, L, R )

c*********************************************************************72
c
C  N  - NUMBER OF NODES TO BE CONNECTED.
C  P  - NODE NUMBER VECTOR (IN OUTPUT, OPTIMAL CONNECTION).
C  SN - START NODE NUMBER.
C  EN - END NODE NUMBER.
C  M  - DISTANCE MATRIX ORDER.
C  D  - DISTANCE MATRIX.
C  L  - SHORTEST CONNECTION LENGTH (OUTPUT).
C  R  - NUMBER OF RUNS.
C  GET LARGE NUMBER ( = N X MAX D(I,J) ).
c
      INTEGER P(N), D(M,M), ID(60), Q(60), SN, EN, R

      LARGE = 0
      DO 20 I = 1, M
        DO 10 J = 1, M
          IF ( D(I,J) .GT. LARGE ) LARGE = D(I,J)
10      CONTINUE
20    CONTINUE
      LARGE = LARGE * N
C  DEFINE NON-EXISTING ARCS BY ASSIGNING
C  THEIR DISTANCES LARGE NEGATIVE VALUES.
      IF ( EN .NE. 0 ) GO TO 40
      DO 30 I = 1, M
        ID(I) = D(I,SN)
        D(I,SN) = -LARGE
        D(SN,SN) = 0
30    CONTINUE
40    IF ( SN .EQ. EN .OR. EN .EQ. 0 ) GO TO 50
      ID(1) = D(EN,SN)
      D(EN,SN) = -LARGE
C  RUN R TRIALS.
50    L = LARGE
      DO 280 IRS = 1, R
C  BUILD TOUR BY SUCCESSIVE INSERTING
C  NOT-YET-INVOLVED NODES.
C  INITIATE TOUR IS CONSIDERED AS 
C  ARC P(1) TO P(1).
        DO 90 JS = 2, N
          MININC = LARGE
C  TRACE ALL NOT-YET-INVOLVED NODES
C  TO CHOOSE THE ONE WITH MINIMUM INCREMENT.
          DO 70 J = JS, N
            JP = P(J)
            JE = JS - 1
C  FOR EACH NOT-YET-INVOLVED NODE TRACE ALREADY
C  BUILT-UP TOUR TO CHOOSE THE MINIMUM INCREMENT ARC.
            DO 60 I = 1, JE
              IP = P(I)
              IP1 = P(I+1)
              IF ( I .EQ. JE ) IP1 = P(1)
              INC = D(IP,JP) + D(JP,IP1) - D(IP,IP1)
              IF ( INC .GE. MININC ) GO TO 60
              J1 = J
              I1 = I
              MININC = INC
60          CONTINUE
70        CONTINUE
C  STRETCH TOUR BY INSERTING THE CHOSEN NODE P(J1)
C  BETWEEN THE NODES P(I1) AND P(I1+1).
80        J1 = J1 - 1
          IF ( J1 .EQ. I1 ) GO TO 90
          IP = P(J1)
          P(J1) = P(J1+1)
          P(J1+1) = IP
          GO TO 80
90      CONTINUE
C  CORRECT TOUR BY 3-OPT METHOD.
C  VARY CONSECUTIVE CHAIN LENGTH K.
        N1 = N - 1
        IF ( N .LT. 3 ) GO TO 210
        DO 200 K = 1, N1
          ICOUNT = 0
C  SHIFT CONSECUTIVE CHAIN
C  THROUGHOUT SEQUENCE OF N NODES.
100       ICOR = 0
          DO 190 J = 1, N
C  CALCULATE CHAIN LENGTH IN FORWARD
C  AND BACKWARD DIRECTION.
            L1 = 0
            LR = 0
            IF ( K .EQ. 1 ) GO TO 120
            I = J
            K1 = 1
110         IF ( I .GT. N ) I = I - N
            IP = P(I)
            IP1 = I + 1
            IF ( IP1 .GT. N ) IP1 = 1
            IP1 = P(IP1)
            L1 = L1 + D(IP,IP1)
            LR = LR + D(IP1,IP)
            I = I + 1
            K1 = K1 + 1
            IF ( K1 .LT. K ) GO TO 110
C  FOR EACH POSITIONED CHAIN (AS IS AND INVERTED)
C  CHECK ALL ARCS IF INSERTION IMPROVES TOUR.
120         MININC = LARGE
            J1 = J + K - 1
            IF ( J1 .GT. N ) J1 = J1 - N
            DO 150 I = 1, N
              IF ( J .LE. J1 .AND. ( I .GE. J .AND. I .LE. J1 ) )
     &          GO TO 150
              IF ( J .GT. J1 .AND. ( I .LE. J1 .OR. I .GE. J ) )
     &          GO TO 150
              IP = P(I)
              JP = P(J)
              JP1 = P(J1)
              IP1 = I + 1
              IF ( IP1 .GT. N ) IP1 = 1
              JE = IP1
              IF ( IP1 .EQ. J ) IP1 = J1 + 1
              IF ( IP1 .GT. N ) IP1 = 1
              IP1 = P(IP1)
              LN = L1
              IR = 0
130           INC = D(IP,JP) + LN + D(JP1,IP1) - D(IP,IP1)
              IF ( INC .GT. MININC .OR. ( INC .EQ. MININC .AND.
     &          ( JE .NE. J .OR. ( JE .EQ. J .AND. IR .EQ. 1 ) ) ) ) 
     &          GO TO 140
              I1 = I
              IR1 = IR
              MININC = INC
140           IF ( IR .EQ. 1 ) GO TO 150
              IR = 1
              LN = LR
              JS = JP
              JP = JP1
              JP1 = JS
              GO TO 130
150         CONTINUE
            I = I1 + 1
            IF ( I .GT. N ) I = 1
            IF ( I .EQ. J .AND. IR1 .EQ. 0 ) GO TO 190
C  REINSERT CHAIN OF LENGTH K STARTING IN J
C  BETWEEN NODES P(I1) AND P(I1+1).
            ICOR = 1
            JS = J
            JE = 0
            IF ( IR1 .EQ. 0 ) GO TO 160
            JS = J1
            JE = -1
160         K1 = 0
170         K1 = K1 + 1
            IF ( K1 .GT. K ) GO TO 190
            I = JS
            JS = JS + JE
            IF ( JS .LT. 1 ) JS = N
180         IP = I + 1
            IF ( IP .GT. N ) IP = 1
            JP = P(I)
            P(I) = P(IP)
            P(IP) = JP
            I = I + 1
            IF ( I .GT. N ) I = 1
            IF ( IP - I1 ) 180, 170, 180
190       CONTINUE
          IF ( ICOR .EQ. 0 ) GO TO 200
          ICOUNT = ICOUNT + 1
          IF ( ICOUNT .LT. N ) GO TO 100
200     CONTINUE
C  ORIENT TOUR WITH SN IN P(1).
210     DO 230 I = 1, N
          IF ( P(1) .EQ. SN ) GO TO 240
          JS = P(1)
          DO 220 J = 1, N1
            P(J) = P(J+1)
220       CONTINUE
          P(N) = JS
230     CONTINUE
C  CALCULATE TOUR LENGTH.
240     L1 = 0
        DO 250 I = 1, N1
          IP = P(I)
          IP1 = P(I+1)
          L1 = L1 + D(IP,IP1)
250     CONTINUE
        IP = P(1)
        IF ( SN .EQ. EN ) L1 = L1 + D(IP1,IP)
C  SAVE SOLUTION, IF BETTER, AND SET NEW INITIATE NODE.
        IF ( L1 .GE. L ) GO TO 270
        L = L1
        DO 260 I = 1,  N
          Q(I) = P(I)
260     CONTINUE
270     J = IRS + 1
        IF ( J .GT. N ) J = J - N
        JS = P(1)
        P(1) = P(J)
        P(J) = JS
280   CONTINUE
C  RESTORE P AND DUMMY DISTANCES.
      DO 290 I = 1, N
        P(I) = Q(I)
290   CONTINUE
      IF ( EN .NE. 0 ) GO TO 310
      DO 300 I = 1, M
        D(I,SN) = ID(I)
300   CONTINUE
310   IF ( SN .EQ. EN .OR. EN .EQ. 0 ) GO TO 320
      D(EN,SN) = ID(1)
320   RETURN
      END
