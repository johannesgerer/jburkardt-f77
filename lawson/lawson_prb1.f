      program main 

c*********************************************************************72
c
cc MAIN is the main program for LAWSON_PRB1.
c
C  DEMONSTRATE ALGORITHMS HFT AND HS1 FOR SOLVING LEAST SQUARES   
C  PROBLEMS AND ALGORITHM COV FOR COMPUTING THE ASSOCIATED COVARIANCE
C  MATRICES. 
C   
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
      integer MDA
      parameter(MDA = 8)
      integer I, IP1, J, JM1, K, L
      integer M, MMN, MN1, MN2, N, NM1, NOISE, NP1, NPJ
      double precision A(MDA,MDA), ANOISE, B(MDA), DUMMY, GEN, H(MDA)
      double precision ONE, SM, SMALL, SRSMSQ, ZERO
      parameter(ONE = 1.0d0, SMALL = 1.0d-4, ZERO = 0.0d0)

  300 format (/1x,8X,'ESTIMATED PARAMETERS,  X=A**(+)*B,',
     * ' COMPUTED BY ''HFT,HS1'''//
     * (9X,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8))  
  305 format (/1x,8X,'RESIDUAL LENGTH =',E12.4)
  310 format (/1x,8X,'COVARIANCE MATRIX (UNSCALED) OF',
     * ' ESTIMATED PARAMETERS'/
     * 9x,'COMPUTED BY ''COV''.'/1X)  
  320 format (9X,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8)     
  330 format (/' PROG1.  THIS PROGRAM DEMONSTRATES THE ALGORITHMS',
     * ' HFT, HS1, AND COV.'//
     * ' CAUTION.. ''PROG1'' DOES NO CHECKING FOR',
     * ' NEAR RANK DEFICIENT MATRICES.'/
     * ' RESULTS IN SUCH CASES MAY BE MEANINGLESS.'/
     * ' SUCH CASES ARE HANDLED BY ''PROG2'' OR ''PROG3''')    
  340 format (/
     * ' THE RELATIVE NOISE LEVEL OF THE GENERATED DATA WILL BE ',g11.3)
  350 format (/////'    M   N'/1X,2I4)
  360 format (/2x,
     * '******  TERMINATING THIS CASE DUE TO',
     * ' A DIVISOR BEING EXACTLY ZERO. ******')

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB1'

      DO 230 NOISE=1,2 
         if(NOISE .eq. 1) then
            ANOISE= ZERO
         else
            ANOISE= SMALL
         endif
         write (*,330)
         write (*,340) ANOISE    
c 
C  INITIALIZE THE DATA GENERATION FUNCTION   
C 
         DUMMY=GEN(-ONE)   
         DO 220 MN1=1,6,5    
            MN2=MN1+2   
            DO 210 M=MN1,MN2   
               DO 200 N=MN1, M
                  NP1=N+1
                  write (*,350) M,N
c
C  GENERATE DATA     
C
                  DO 10 I=1,M   
                     DO 10 J=1,N  
   10                   A(I,J)=GEN(ANOISE)   
                  DO 20 I=1,M   
   20                B(I)=GEN(ANOISE)  
                  IF(M .LT. N) GO TO 180     
C   
C     ******  BEGIN ALGORITHM HFT  ******   
C    ..     
                  DO 30 J=1,N   
                     CALL H12 (1,J,J+1,M,A(1,J),1,H(J),
     *                    A(1,min(J+1,N)),1,MDA,N-J)
   30             continue
C    ..     
C     THE ALGORITHM 'HFT' IS COMPLETED. 
C   
C     ******  BEGIN ALGORITHM HS1  ******   
C     APPLY THE TRANSFORMATIONS  Q(N)...Q(1)=Q TO B 
C     REPLACING THE PREVIOUS CONTENTS OF THE ARRAY, B .     
C    ..     
                  DO 40 J=1,N   
   40                CALL H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,1,1)
C     SOLVE THE TRIANGULAR SYSTEM FOR THE SOLUTION X.   
C     STORE X IN THE ARRAY, B .     
C    ..     
                  DO 80 K=1,N   
                     I=NP1-K   
                     SM= ZERO
                     IF (I.EQ.N) GO TO 60  
                     IP1=I+1   
                     DO 50 J=IP1,N    
   50                   SM=SM+A(I,J)*B(J)
   60                IF (A(I,I) .eq. ZERO) then
                         write (*,360) 
                         GO TO 180 
                     endif
                     B(I)=(B(I)-SM)/A(I,I) 
   80             continue
c
C  COMPUTE LENGTH OF RESIDUAL VECTOR.   
C
                  SRSMSQ= ZERO
                  IF (N.EQ.M) GO TO 100  
                  MMN=M-N
                  DO 90 J=1,MMN 
                     NPJ=N+J   
   90                SRSMSQ=SRSMSQ+B(NPJ)**2   
                  SRSMSQ=SQRT(SRSMSQ)
C     ******  BEGIN ALGORITHM  COV  ******  
C     COMPUTE UNSCALED COVARIANCE MATRIX   ((A**T)*A)**(-1) 
C
  100             DO 110 J=1,N  
  110                A(J,J)= ONE/A(J,J)  
                  IF (N.EQ.1) GO TO 140  
                  NM1=N-1
                  DO 130 I=1,NM1
                     IP1=I+1   
                     DO 130 J=IP1,N   
                        JM1=J-1  
                        SM= ZERO
                        DO 120 L=I,JM1  
  120                      SM=SM+A(I,L)*DBLE(A(L,J))   
  130                   A(I,J)=-SM*A(J,J)
C
C  THE UPPER TRIANGLE OF A HAS BEEN INVERTED UPON ITSELF.  
c
  140             DO 160 I=1,N  
                     DO 160 J=I,N     
                        SM= ZERO 
                        DO 150 L=J,N
  150                      SM=SM+A(I,L)*DBLE(A(J,L))   
  160                   A(I,J)=SM
C
C  THE UPPER TRIANGULAR PART OF THE SYMMETRIC MATRIX (A**T*A)**(-1) HAS   
C  REPLACED THE UPPER TRIANGULAR PART OF THE A ARRAY.  
c
                  write (*,300) (I,B(I),I=1,N)   
                  write (*,305) SRSMSQ   
                  write (*,310)  
                  DO 170 I=1,N  
  170                write (*,320) (I,J,A(I,J),J=I,N)  
  180             CONTINUE   
  200          continue
  210       continue
  220    continue
  230 continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB1:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      STOP  
      END   
