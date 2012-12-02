      PROGRAM MAIN

c*********************************************************************72
c
CC MAIN is the main program for KMEDIAN.
C
C             LAGRANGIAN RELAXATION METHOD FOR
C             K-MEDIAN PROBLEM  (WITH WEIGHT)
C          ALGORITHM OF CORNUEJOLS,FISHER,NEMHAUSER
C 
C
      INTEGER UNIT
      COMMON/IOLIST/IDISTA(50,50),IDEND(50)
      COMMON/LL/ISHODI(50,50),NO,KOPEN
      COMMON/SRT/JSRT(50,50)
      COMMON/PERIP/UNIT
      COMMON/SOLMU/NVSQ,NV(50),W,ITERA,NREAL,BESTW
      COMMON/HEURIS/NRWARD,JWARHS(50),IS(50)
      COMMON/LAGRAN/U(50),ISW
C
      OPEN(31,FILE='KMEDI.DAT',STATUS='NEW')
C     
C     CALL I/O SUBROUTINE
C
10    CONTINUE
C
      CALL KWRITER
C
C
5892  DO 5916 I = 1,NO
              DO 5917 J = 1,NO
                 ISHODI(I,J) = IDISTA(I,J)*IDEND(I)
5917          CONTINUE
5916  CONTINUE
5819  ITERA = 0
      ISW = 0
      CALL GREEDY
      BESTW = -1000000.
      ITER = 0
      AMBDK = 2.
      DO 13 I = 1,NO
                DO 17 J = 1,NO
                      IS(J)=ISHODI(I,J)
                      IDEND(J)=J
17              CONTINUE
      CALL SORT(IS,IDEND,NO)
                DO 18 J= 1,NO
                      JSRT(I,J)=IDEND(J)
18              CONTINUE
      U(I)=-ISHODI(I,IDEND(1))    
13    CONTINUE
C
C       SOLUTION OF FIXED LAGRANGIAN MULTIPLIER
C
14    CALL SOL
      IF(ISW .NE. 0) GO TO 1000
      IF(W .LE. BESTW) GO TO 45
      BESTW = W
      ITER = 0
      GO TO 15
45    IF(ITER .LT. 8) GO TO 15
      AMBDK = AMBDK/2
      ITER = 0
C
C     COMPUTATION OF A NEW MULTIPLIER VECTOR  FOR SUBGRADIENT ITERATION
C
15    TK = (AMBDK*(NRWARD-W))/NVSQ
      ITER = ITER + 1
      DO 16 I =1,NO
16    U(I) = U(I) + TK*NV(I)
      WRITE(31,864) NRWARD,BESTW
864   FORMAT(' BEST HEURISTIC SOLUTION',1X,I8,
     X '     BEST LAGRANGIAN BOUND',1X,F12.2)
      GO TO 14
1000  CONTINUE
      WRITE(*,*) ' YOUR OUTPUT IS IN FILE KMEDI.DAT '
      STOP
      END
      SUBROUTINE SOL

c*********************************************************************72
C
CC SOL takes one step of the iteration.
C
      COMMON/IOLIST/IDISTA(50,50),IDEND(50)
      COMMON/LL/ISHODI(50,50),NO,KOPEN
      COMMON/SRT/JSRT(50,50)
      COMMON/PERIP/UNIT
      COMMON/SOLMU/NVSQ,NV(50),W,ITERA,NREAL,BESTW
      COMMON/LAGRAN/U(50),ISW
      COMMON/HEURIS/NRWARD,JWARHS(50),IS(50)
      DIMENSION ROSUM(50),ICJMAX(50),ROWJ(50),KS(50)
      DO 10 J = 1,NO
                ROSUM(J) = 0
10    CONTINUE
      DO 12 I = 1,NO
      ITR=1
11    NEXTJ=JSRT(I,ITR)
      SUMIND = ISHODI(I,NEXTJ) + U(I)
      IF(SUMIND.GE.0) GOTO 12
      ROSUM(NEXTJ) = ROSUM(NEXTJ) + SUMIND
      ITR=ITR+1
      IF(ITR.LE.NO) GOTO 11
12    CONTINUE
      DO 771 J=1,NO
771   ROWJ(J) = ROSUM(J)
      DO 14 K = 1,KOPEN
            COMPAR = 1000000.
                 DO 13 J = 1,NO
                    IF(COMPAR .LE. ROWJ(J)) GO TO 13
                    COMPAR = ROWJ(J)
                    INDXJ = J
13               CONTINUE
            ICJMAX(K) =INDXJ
            ROWJ(INDXJ) = 1000000.
14    CONTINUE
      WW = 0.
      DO 23 L = 1,KOPEN
            IODX = ICJMAX(L)
            WW = WW + ROSUM(IODX)
23    CONTINUE
      DO 30 I = 1,NO
            WW = WW - U(I)
30    CONTINUE
      W = WW
      IF(W .GT. BESTW-0.001 .AND. W .LT. BESTW+0.001) GO TO 72
      INDX = 0
      GO TO 73
72    IF(KOPEN .EQ. NO-1) GO TO 74
      INDX = INDX + 1
73    IF(INDX .LT. 20) GO TO 74
      WRITE(31,374)
374   FORMAT(' THE LARGRANGIAN RELAXATION CONVERGES TO A SPECIFIC
     1 VALUE ')
      ISW = 3
      GO TO 27
74    IF(KOPEN .EQ. 0) GO TO 63
      NREAL = 0
      DO 45 I = 1,NO
            JJMAX = 1000000
                  DO 46 JJ = 1,KOPEN
                     JK = ICJMAX(JJ)
                     IF(ISHODI(I,JK) .LT. JJMAX)JJMAX =ISHODI(I,JK)
46                CONTINUE
            NREAL = NREAL + JJMAX
45    CONTINUE
      IF(NRWARD .LE. NREAL) GO TO 63
      NRWARD = NREAL
      DO 3431 KKK = 1,KOPEN
3431  JWARHS(KKK) = ICJMAX(KKK)
63    IF(BESTW + 0.9 .GE. NRWARD) GO TO 400
      GO TO 401
400   CONTINUE
      ISW = 1
      GO TO 27
401   DO 25 I = 1,NO
          NV(I) = -1
          IF(KOPEN .EQ. 0) GO TO 25
               DO 24 J = 1,KOPEN
                  IOIDX = ICJMAX(J)
                  IF(ISHODI(I,IOIDX)+U(I) .LE. 0.) NV(I) =NV(I) + 1
24             CONTINUE
25    CONTINUE
      NVSQ = NV(1)**2
      DO 26 JJ = 2,NO
               NVSQ = NVSQ + NV(JJ)**2
26    CONTINUE
      ITERA = ITERA + 1
      IF(NVSQ .EQ. 0) GO TO 500
      GO TO 501
500   continue  
      ISW = 2
      GO TO 27
501   IF(ITERA .LE. 300) GO TO 29
      ISW = 4
      IF(W .GT. BESTW) BESTW = W
      WRITE(31,*) 'THE COST OF THE BEST SOLUTION FOUND IS:',NRWARD
      GO TO 289
27    IF(W .GT. BESTW) BESTW = W  
289   CONTINUE  
      IF(ISW .EQ. 1)GOTO 236
      IF(ISW .EQ. 2)GOTO 236
      WRITE(31,235)NRWARD
235   FORMAT(/,' BEST FEASIBLE SOLUTION',  I12)
      GOTO 511
236   WRITE(31,237)NRWARD
237   FORMAT(/,'COST OF OPTIMAL SOLUTION IS : ', I12)
511   CONTINUE
      WRITE(31,*)
      WRITE(31,104)(JWARHS(J),J=1,KOPEN)
104   FORMAT(' THE STOCKS IN THE INDEX FUND ARE:'/20(1X,I4))
      WRITE(31,*)     
      WRITE(31,*)'FUND STOCK  TRACKS  TARGET STOCKS'
      DO 129 I = 1,NO
         MCOMP = 10000000
         DO 229 J = 1,KOPEN
             IF(MCOMP.LE.ISHODI(I,JWARHS(J)))GO TO 229
             MCOMP=ISHODI(I,JWARHS(J))
             KS(I) = JWARHS(J)         
229      CONTINUE    
      WRITE(*,*)
      WRITE(31,329)KS(I),I
329   FORMAT(I3,'            ----->    ',I5)
129   CONTINUE
      IF(ITERA.LT.300) GOTO 512
      WRITE(31,*)' AFTER 300 ITERATIONS LAGRANGIAN BOUND IS:',BESTW
      GOTO 29
512   WRITE(31,*)' NUMBER OF ITERATIONS:',ITERA
29    RETURN
      END
      SUBROUTINE GREEDY

c*********************************************************************72
C
CC GREEDY finds a nonoptimal solution using a greedy algorithm.
C
      COMMON/IOLIST/IDISTA(50,50),IDEND(50)
      COMMON/LL/ISHODI(50,50),NO,KOPEN
      COMMON/PERIP/UNIT
      COMMON/HEURIS/NRWARD,JWARHS(50),IS(50)
      DIMENSION IUU(50),ICOLM(50)
      ICOLTO = 0
      NWAR = 0
      DO 5 I = 1,NO
5     IS(I) = 0
      IZVA = 0
      DO 30 I = 1,NO
      IUMAX = -10**5
      DO 20 J = 1,NO
      IF(ISHODI(I,J) .LT.IUMAX) GO TO 20
      IUMAX =ISHODI(I,J)
20    CONTINUE
      IUU(I) = IUMAX
      IZVA = IZVA + IUU(I)
30    CONTINUE
100   DO 150 JS = 1,NO
150   ICOLM(JS) =0
      DO 200 J = 1,NO
      IF(IS(J) .EQ. 1) GO TO 200
      ICOLSM = 0
      DO 250 I =1,NO
      IF(IUU(I) - ISHODI(I,J) .LT. 0) GO TO 250
      ICOLSM=ICOLSM+IUU(I)-ISHODI(I,J)
250   CONTINUE
      ICOLM(J) = ICOLSM
200   CONTINUE
      KYU = -100000
      DO 260 JJS = 1,NO
      IF(IS(JJS) .EQ. 1) GO TO 260
      IF(ICOLM(JJS) .LT. KYU) GO TO 260
      KYU = ICOLM(JJS)
      MAXIDX = JJS
260   CONTINUE
      ICOLTO = ICOLTO + KYU
      IS(MAXIDX) = 1
      NWAR = NWAR + 1
      JWARHS(NWAR) = MAXIDX    
      IF(NWAR .EQ. KOPEN) GO TO 400
      DO 300 I = 1,NO
      IF(IUU(I) - ISHODI(I,MAXIDX) .LE. 0) GO TO 300
      IUU(I) = ISHODI(I,MAXIDX)
300   CONTINUE
      GO TO 100
400   NRWARD = IZVA - ICOLTO
      RETURN
      END
      SUBROUTINE KWRITER

c*********************************************************************72
c
CC KWRITER gets input from the user.
C
      COMMON/LL/ISHODI(50,50),NO,KOPEN
      COMMON/IOLIST/IDISTA(50,50),IDEND(50)
      CHARACTER*64 FINNAME
      CHARACTER*64 FOUTNAME
C
      WRITE(*,*) '*************************************************'
      WRITE(*,*) '*                                               *'
      WRITE(*,*) '*       LAGRANGIAN RELAXATION METHOD FOR        *'
      WRITE(*,*) '*                                               *'
      WRITE(*,*) '*       THE K-MEDIAN PROBLEM WITH WEIGHTS.      *'
      WRITE(*,*) '*                                               *'
      WRITE(*,*) '*                 ALGORITHM OF                  *'
      WRITE(*,*) '*                                               *'
      WRITE(*,*) '*       CORNUEJOLS, FISHER AND NEMHAUSER.       *'
      WRITE(*,*) '*                                               *'
      WRITE(*,*) '*                                         1977  *'
      WRITE(*,*) '*************************************************'
      WRITE(*,*)
      WRITE(*,*) 'DO YOU HAVE A SAVED FILE? ( 0 => YES, >0 => NO):'
      READ(*,*) KSAVED
      IF ( KSAVED .EQ. 0 ) GOTO 1010       
      WRITE(*,*)
      WRITE(*,*) 'ENTER NO OF CITIES ( NO =< 50 ): '
      READ(*,*) NO
      WRITE(*,*)
      WRITE(*,*) 'ENTER NO OF OPEN FACILITIES ( KOPEN =< 50 ):'
      READ(*,*) KOPEN
C
C     INITIALIZE DISTANCE MATRIX, IDISTA(I,J) = 9999999.
C
      DO 1000 I1 = 1,NO
         DO 1001 I2 = 1,NO
            IDISTA(I1,I2) = 99999
1001     CONTINUE
1000  CONTINUE
C
      WRITE(*,*) 'DISTANCE MATRIX IS CURRENTLY FILLED WITH DEFAULT '
      WRITE(*,*) 'VALUES: 99,999. ENTER ONLY THE ONES YOU WISH TO '
      WRITE(*,*) 'CHANGE. ENTER THREE PARAMETERS [ I, J, C(I,J) ]:'
      WRITE(*,*) 'ENTER 0, 0, 0, TO STOP.'      
100   CONTINUE        
      READ(*,*) IR,IC,CIJ
      IF( IR . LE. 0 ) GOTO 101
      IDISTA(IR,IC) = CIJ
      IDISTA(IC,IR) = CIJ
      GOTO 100      
101   CONTINUE
C
      WRITE(*,*) 'ENTER THE WEIGHT OF STOCK:'
      DO 1003 I1 = 1,NO
         WRITE(*,*) 'STOCK (',I1,')'
         READ(*,*) IDEND(I1)
1003  CONTINUE
      WRITE(*,*) ' WISH TO SAVE INPUT FILE? ( 0 YES, >0 NO):'
      READ(*,*) KSAVE 
      IF(KSAVE .GT. 0 ) go to 1504
      WRITE(*,*) ' ENTER THE FILENAME:'
      READ(*,1011) FOUTNAME
      OPEN(22,FILE=FOUTNAME,STATUS='NEW')
      WRITE(22,1012)NO,KOPEN
      DO 1501 II=1,NO
        WRITE(22,1013) IDEND(II)
1501  CONTINUE
      DO 1502 I1 = 1,NO
        DO 1503 I2 = 1,NO  
          WRITE(22,1013) IDISTA(I1,I2)
1503  CONTINUE
1502  CONTINUE
1504  CONTINUE
      GOTO 2000
C
1010  WRITE(*,*) ' INPUT FILE NAME:'
      READ(*,1011) FINNAME
C     
C     OPEN UNIT 21 FOR INPUT SAVED DATA ...
C
      OPEN(21,FILE=FINNAME)
      READ(21,1012)NO,KOPEN
      DO 1020 II=1,NO 
        READ(21,1013) IDEND(II)
1020  CONTINUE
      DO 1030 I1=1,NO
        DO 1040 I2=1,NO
          READ(21,1013) IDISTA(I1,I2)
1040    CONTINUE  
1030  CONTINUE
C
C     FORMATS...
C
1011  FORMAT(A)
1012  FORMAT(1X,2I4)
1013  FORMAT(1X,I5)
2000  CONTINUE            
      RETURN
      END
      SUBROUTINE  SORT(RR,PT,N)

c*********************************************************************72
c
CC SORT sorts the closest fund stock for each target stock.
C
      INTEGER RR,PT
      DIMENSION RR(50),PT(50)
      I=1
10    M=I+I+1
      I=I+1
      IF(N.GE.I) GOTO 10
14    M=M/2
      IF(M.EQ.0) RETURN
18    K=N-M
      J=1
20    IF(J.GT.K) GOTO 14
      I=J
24    IF(I.GE.1) GOTO 28
26    J=J+1
      GOTO 20
28    MM=M+I
      IF(RR(MM).GE.RR(I)) GOTO 26
30    IT4=PT(I)
      PT(I)=PT(MM)
      PT(MM)=IT4
      IT3=RR(I)
      RR(I)=RR(MM)
      RR(MM)=IT3
      I=I-M 
      GOTO 24
      END
        
