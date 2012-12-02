      program main

c*********************************************************************72
c
cc MAIN is the main program for LAWSON_PRB2.
c
c  Discussion:
c
C    Demonstrate algorithm HFTI for solving least squares problems
c    and algorithm COV for computing the associated unscaled 
c    covariance matrix.
c
c  Modified:
c
c    22 October 2008
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
      parameter ( MDA = 8 )

      double precision A(MDA,MDA)
      double precision ANOISE
      double precision ANORM
      double precision B(MDA)
      double precision DUMMY
      double precision FIVE00
      double precision G(MDA)
      double precision GEN
      double precision H(MDA)
      double precision HALF
      integer I
      integer II
      integer IP(MDA)
      integer IP1
      integer J
      integer JM1
      integer K
      integer KM1
      integer KP1
      integer KRANK
      integer L
      integer M
      integer MN1
      integer MN2
      integer N
      integer NM1
      integer NOISE
      double precision ONE
      double precision SM
      double precision SMALL
      double precision SRSMSQ
      double precision TAU
      double precision TMP

      parameter(FIVE00 = 500.0d0, HALF = 0.5d0)
      parameter(ONE = 1.0d0, SMALL = 1.0d-4 )
  
  190 format (/1x,8X,'RESIDUAL LENGTH = ',E12.4)
  200 format (/1x,8X,
     * 'ESTIMATED PARAMETERS,  X=A**(+)*B, COMPUTED BY ''HFTI'''//
     * (9X,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8,I6,E16.8))  
  210 format (/1x,8X,
     * 'COVARIANCE MATRIX (UNSCALED) OF ESTIMATED PARAMETERS'/
     * 9x,'COMPUTED BY ''COV''.'/1X)  
  220 format (9X,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8,2I3,E16.8)     
  250 format (/////'    M   N'/1X,2I4)
  260 format (/1x,8X,'PSEUDORANK = ',I4) 

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB2'
      write ( *, '(a)' ) '  Demonstrate the algorithms HFTI and COV.'

      DO 180 NOISE = 1, 2  

        ANORM = 500.0D+00

        if ( noise .eq. 1 ) then
          ANOISE = 0.0D+00
          TAU = HALF
        else
          ANOISE = SMALL
          TAU = ANORM * ANOISE * 10.0D+00
        end if
c
C  INITIALIZE THE DATA GENERATION FUNCTION   
C     
        DUMMY=GEN(-ONE)

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) 
     &    '  Use a relative noise level of        ', anoise
        write ( *, '(a,g14.6)' ) 
     &    '  The matrix norm is approximately     ', anorm
        write ( *, '(a,g14.6)' ) 
     &    '  The absolute pseudorank tolerance is ', tau
   
        DO 180 MN1 = 1, 6, 5  
          MN2 = MN1 + 2 
          DO 180 M = MN1, MN2  
            DO 180 N = MN1, MN2  

              write ( *, '(a)' ) ' '
              write ( *, '(a,i6)' ) '  M = ', m
              write ( *, '(a,i6)' ) '  N = ', n
c
C  GENERATE DATA     
C     
              do i = 1, m   
                do j = 1, n   
                  a(i,j) = gen ( anoise )
                end do
              end do

              do i = 1, m   
                b(i) = gen ( anoise )  
              end do
C   
C  CALL HFTI.   
C   
              CALL HFTI(A,MDA,M,N,B,1,1,TAU,KRANK,SRSMSQ,H,G,IP)

              write (*,260) KRANK   
              write (*,200) (I,B(I),I=1,N)  
              write (*,190) SRSMSQ  

              IF (KRANK.LT.N) GO TO 180 
c
C  ALGORITHM COV BEGINS HERE.
C  
              DO J=1,N   
                A(J,J)= ONE/A(J,J)  
              end do

                      IF (N.EQ.1) GO TO 70  
                      NM1=N-1   
                          DO 60 I=1,NM1 
                          IP1=I+1   
                              DO 60 J=IP1,N     
                              JM1=J-1   
                              SM= 0.0D+00
                                  DO 50 L=I,JM1 
   50                             SM=SM+A(I,L)*A(L,J)
   60                         A(I,J)=-SM*A(J,J) 
C
C  THE UPPER TRIANGLE OF A HAS BEEN INVERTED UPON ITSELF.
c
   70                     DO 90 I=1,N   
                              DO 90 J=I,N   
                              SM= 0.0D+00   
                                  DO 80 L=J,N   
   80                             SM=SM+A(I,L)*DBLE(A(J,L)) 
   90                         A(I,J)=SM 

                      IF (N.LT.2) GO TO 160     
                          DO 150 II=2,N 
                          I=N+1-II  
                          IF (IP(I).EQ.I) GO TO 150 
                          K=IP(I)   
                          TMP=A(I,I)
                          A(I,I)=A(K,K) 
                          A(K,K)=TMP
                          IF (I.EQ.1) GO TO 110 
                              DO 100 L=2,I  
                              TMP=A(L-1,I)  
                              A(L-1,I)=A(L-1,K) 
  100                         A(L-1,K)=TMP  
  110                     IP1=I+1   
                          KM1=K-1   
                          IF (IP1.GT.KM1) GO TO 130 
                              DO 120 L=IP1,KM1  
                              TMP=A(I,L)
                              A(I,L)=A(L,K)     
  120                         A(L,K)=TMP
  130                     IF (K.EQ.N) GO TO 150 
                          KP1=K+1   
                              DO 140 L=KP1,N    
                              TMP=A(I,L)
                              A(I,L)=A(K,L)     
  140                         A(K,L)=TMP
  150                     CONTINUE  
  160                 CONTINUE  
C
C  COVARIANCE HAS BEEN COMPUTED AND REPERMUTED.  THE UPPER TRIANGULAR PART OF THE  
C  SYMMETRIC MATRIX (A**T*A)**(-1) HAS REPLACED THE UPPER TRIANGULAR PART OF     
C  THE A ARRAY.
c
                      write (*,210) 
                          DO I=1,N  
                          write (*,220) (I,J,A(I,J),J=I,N)  
                          end do

  180                 CONTINUE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAWSON_PRB2:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
